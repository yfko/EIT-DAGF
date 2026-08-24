function [imf, res, info] = dagaf2(S, nIMF, chi, eType1, eType2, TH1, TH2, alpha, maskCap)
% DAGAF2  Data-adaptive Gaussian average filtering (silent, corrected).
%
% Faithful re-implementation of:
%   Y.-D. Lin, Y. K. Tan, B. Tian, "A novel approach for decomposition of
%   biomedical signals in different applications based on data-adaptive
%   Gaussian average filtering", Biomed. Signal Process. Control 71:103104, 2022.
%
% Differences from the authors' released dagaf.m (all deliberate):
%   1. alpha defaults to 3 (paper Eq.6), not 4.0728.
%   2. Energy ratio follows Eq.(19) exactly: 10*log10(||s||^2/||r||^2).
%      The released code used 10*log10(||s||/||r||), i.e. half the dB value.
%   3. eType1/eType2 are honoured. The released code's validation test
%      (a ~= 'p' || a ~= 'c' || ...) is true for every input, so the
%      extension type was always silently forced to 'd'.
%   4. No figures are drawn (this is called thousands of times in a sweep).
%   5. gausswin() is inlined, so no Signal Processing Toolbox is required.
%   6. maskCap lets you bound the Gaussian window length -- needed when the
%      "signal" is an EIT measurement vector whose samples are not a
%      uniformly-sampled time series (see eit_dagaf_pipeline.m).
%
% Inputs
%   S        signal (vector)
%   nIMF     number of IMFs to extract                    (default 2)
%   chi      window-length scalar, Eq.(17), in [1.1, 3]   (default 2)
%   eType1   left  boundary extension 'p'|'c'|'r'|'d'     (default 'd')
%   eType2   right boundary extension 'p'|'c'|'r'|'d'     (default 'd')
%   TH1      stop when energy ratio (dB) exceeds this     (default 20)
%   TH2      stop when relative change falls below this   (default 1e-3)
%   alpha    Gaussian window shape, Eq.(5)-(6)            (default 3)
%   maskCap  hard upper bound on M (Inf = paper behaviour)(default Inf)
%
% Outputs
%   imf   K x N, IMF k in row k, ordered fast -> slow
%   res   1 x N, S - sum(imf,1)   (lowest-frequency residual)
%   info  struct with the M and Ne actually used at each stage

if nargin < 2 || isempty(nIMF),   nIMF   = 2;     end
if nargin < 3 || isempty(chi),    chi    = 2;     end
if nargin < 4 || isempty(eType1), eType1 = 'd';   end
if nargin < 5 || isempty(eType2), eType2 = 'd';   end
if nargin < 6 || isempty(TH1),    TH1    = 20;    end
if nargin < 7 || isempty(TH2),    TH2    = 1e-3;  end
if nargin < 8 || isempty(alpha),  alpha  = 3;     end
if nargin < 9 || isempty(maskCap),maskCap= Inf;   end

if ~any(eType1 == 'pcrd'), eType1 = 'd'; end
if ~any(eType2 == 'pcrd'), eType2 = 'd'; end
if chi > 3 || chi < 1.1,   chi = 2;   end

L  = numel(S);
S  = reshape(double(S), 1, L);
imf = zeros(0, L);
S1  = S;                       % current residual to be decomposed
info = struct('M', [], 'Ne', [], 'stop', 'nIMF reached');

for k = 1:nIMF
    % --- stop criterion 1: energy ratio, Eq.(19) -------------------------
    energyRatio = 10*log10( (norm(S,2)/norm(S1,2))^2 );
    if energyRatio > TH1
        info.stop = sprintf('energy ratio %.1f dB > TH1 at stage %d', energyRatio, k);
        break
    end

    % --- adaptive window length, Eq.(17) ---------------------------------
    Ne = numel(local_extrema(S1));
    if Ne <= 2
        info.stop = sprintf('too few extrema (%d) at stage %d', Ne, k);
        break
    end
    M = 2*floor(chi*L/Ne);
    M = max(1, min([M, floor((L-1)/2), maskCap]));   % Eq.(18) guard + user cap

    % --- normalised truncated Gaussian window, Eq.(5) & (8) --------------
    m = -M:M;
    W = exp(-0.5*(alpha*m/M).^2);        % identical to gausswin(2*M+1, alpha)
    W = W / sum(W);

    % --- boundary extension, Eq.(10)-(13) --------------------------------
    Se = local_extend(S1, M, eType1, eType2);

    % --- instantaneous mean, Eq.(14). W is symmetric so conv == xcorr ----
    ave = conv(Se, W, 'valid');          % length L

    % --- prototype IMF, Eq.(15) ------------------------------------------
    S2 = S1 - ave;

    % --- stop criterion 2: convergence check -----------------------------
    rTol = (norm(S1 - S2, 2) / norm(S1, 2))^2;
    if rTol < TH2
        info.stop = sprintf('converged (rTol %.2e < TH2) at stage %d', rTol, k);
        break
    end

    imf(end+1, :) = S2;                  %#ok<AGROW>
    S1 = ave;                            % next stage works on the smooth part
    info.M(end+1)  = M;                  %#ok<AGROW>
    info.Ne(end+1) = Ne;                 %#ok<AGROW>
end

res = S - sum(imf, 1);
end
% =========================================================================
function Se = local_extend(S, M, e1, e2)
% Eq.(10) constant, (11) periodical, (12) reflection, (13) double-symmetric.
mu = mean(S);
switch e1
    case 'p', left = S(end-M:end-1);
    case 'r', left = fliplr(S(2:M+1));
    case 'd', left = 2*mu - fliplr(S(2:M+1));
    otherwise, left = S(1)*ones(1, M);            % 'c'
end
switch e2
    case 'p', right = S(2:M+1);
    case 'r', right = fliplr(S(end-M:end-1));
    case 'd', right = 2*mu - fliplr(S(end-M:end-1));
    otherwise, right = S(end)*ones(1, M);         % 'c'
end
Se = [left S right];
end
% =========================================================================
function idx = local_extrema(S)
% Indices of local maxima and minima (second-derivative test).
Mx = find(diff(diff(S) > 0) < 0);
j  = find(S(Mx+1) > S(Mx));  Mx(j) = Mx(j)+1;
X  = -S;
Mn = find(diff(diff(X) > 0) < 0);
j  = find(X(Mn+1) > X(Mn));  Mn(j) = Mn(j)+1;
idx = sort([Mx Mn]);
end
