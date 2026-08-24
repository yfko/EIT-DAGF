function [comp, info] = dagaf_axis(X, axis_, nIMF, chi, eType, TH1, TH2, alpha, cap)
%DAGAF_AXIS  Run DAGAF along a chosen axis of the 13x16 measurement matrix.
%
% The 208 measurements of one adjacent-drive frame are 16 injections x 13
% measurements. Filtering along the concatenated 208-point vector crosses 15
% internal discontinuities, one at every injection boundary.
%
%   axis_ = 'chan'  : the 208-point concatenation (what the paper's EMD does)
%   axis_ = 'inj'   : each of the 13 measurement offsets treated as a
%                     16-sample signal across injection angle
%
% The 'inj' axis is physically continuous: row i tracks how one measurement
% offset responds as the drive pair rotates around the tank, and it is exactly
% periodic with period 16, so periodic extension is correct rather than merely
% plausible.
%
% Returns comp: (K+1) x 208, rows summing exactly to X (last row = residual).

if nargin < 3 || isempty(nIMF),  nIMF  = 2;    end
if nargin < 4 || isempty(chi),   chi   = 2;    end
if nargin < 5 || isempty(eType), eType = 'p';  end
if nargin < 6 || isempty(TH1),   TH1   = 40;   end
if nargin < 7 || isempty(TH2),   TH2   = 1e-8; end
if nargin < 8 || isempty(alpha), alpha = 3;    end
if nargin < 9 || isempty(cap),   cap   = Inf;  end

X = reshape(double(X), 1, []);

switch axis_
    case 'chan'
        [imf, res, in] = dagaf2(X, nIMF, chi, eType, eType, TH1, TH2, alpha, cap);
        comp = [imf; res];
        info = in;

    case 'inj'
        A = reshape(X, 13, 16);          % column j = injection j
        K = nIMF;
        C = zeros(K+1, 13, 16);
        info.M = []; info.Ne = [];
        for i = 1:13                     % each measurement offset, across angle
            [imf, res, in] = dagaf2(A(i,:), nIMF, chi, eType, eType, ...
                                    TH1, TH2, alpha, cap);
            k = size(imf,1);
            C(1:k, i, :) = permute(imf, [1 3 2]);
            C(K+1, i, :) = res;          % rows k+1..K stay zero if it stopped early
            info.M(end+1,1:numel(in.M))   = in.M;   %#ok<AGROW>
            info.Ne(end+1,1:numel(in.Ne)) = in.Ne;  %#ok<AGROW>
        end
        comp = reshape(C, K+1, 208);

    otherwise
        error('axis must be ''chan'' or ''inj''');
end
end
