%% EIT_REVISION_R1
% Computations required by the Round-1 revision roadmap (REVIEW_R1).
%
%   R1(a) independent component-selection criterion, not PE/SD/RNG
%   R2    full PE/SD/RNG triple for DAGAF at window caps 3, 4 and 6
%   R3    positive control: plain fixed linear high-pass filters, both axes
%   R4    variability of PE, SD and RNG across forward-mesh density
%   R6    image of the removed component alone, per method
%   R12   sigma_heart 1.5 vs 2.0 against the published triple
%
% One GREIT model is built and reused throughout; the forward mesh is rebuilt
% only where the experiment requires it (R4).
%
% Writes revision_r1_results.txt and revision_r1_results.mat.

clc; clear; close all;
diary off;
if exist('revision_r1_results.txt','file'), delete('revision_r1_results.txt'); end
diary revision_r1_results.txt

t_all = tic;

%% ---- parameters (identical to eit_compare_table1) ----------------------
P.n_elecs   = 16;
P.r_lung    = 0.6;
P.r_heart   = 0.3;
P.sig_lung  = 0.5;
P.sig_heart = 2.0;
P.maxsz     = 0.2;

P.pca_keep     = 1:12;
P.pca_lung_key = 4;
P.emd_nIMF     = 6;
P.emd_lung_key = 1;
P.dagaf = struct('nIMF',6,'chi',2,'alpha',3,'eType','p', ...
                 'maskCap',6,'TH1',40,'TH2',1e-8);

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here,'rcada_emd'));
addpath(fullfile(here,'github_release','src'));

stim = mk_stim_patterns(P.n_elecs,1,[0,1],[0,1],{'no_meas_current'},1);

tgt_heart = [0; 0.1; 1; P.r_heart];
ridx      = [2 4 5];        % eval_GREIT returns AR PE RES SD RNG

%% ---- GREIT model (built once) -----------------------------------------
fprintf('Building GREIT model ...\n');
opt   = struct('imgsz',[64 64],'distr',0,'Nsim',1000, ...
               'target_size',0.02,'noise_figure',0.5);
fmdlG = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1]);
fmdlG.stimulation = stim;
fmdlG = mdl_normalize(fmdlG,1);
imdl  = mk_GREIT_model(mk_image(fmdlG,1),0.2,[],opt);

ctrT = interp_mesh(imdl.rec_model);
ctrR = (ctrT(1:2:end,:) + ctrT(2:2:end,:))/2;
roi_h = sum((ctrR(:,1:2) - [0  0.1]).^2,2) < P.r_heart^2;
roi_l = sum((ctrR(:,1:2) - [ 1 0.1]).^2,2) < P.r_lung^2 | ...
        sum((ctrR(:,1:2) - [-1 0.1]).^2,2) < P.r_lung^2;
roi_b = ~roi_h & ~roi_l & sum(ctrR(:,1:2).^2,2) < 3.8;   % background inside tank
fprintf('  pixels: %d  (heart %d, lung %d, background %d)\n', ...
        numel(roi_h), sum(roi_h), sum(roi_l), sum(roi_b));

%% ---- baseline forward model -------------------------------------------
[fmdl,vh,vi,X] = build_fwd(P, P.maxsz, P.sig_heart, stim);
fprintf('  %d measurements per frame\n', numel(X));

D = decompose_all(X, P);

%% =======================================================================
%  R1(a)  independent selection criterion
%  The lung-only and heart-only forward solutions give the measurement-domain
%  signature of each structure. Correlating every component against those two
%  references identifies which structure the component carries, using ground
%  truth available on the phantom and touching neither the reconstruction nor
%  any GREIT figure of merit. It is therefore independent of the quantities
%  later reported.
%  =======================================================================
fprintf('\n');
banner('R1(a)  component selection by an independent criterion');
[~,~,~,X_L] = build_fwd(P, P.maxsz, P.sig_heart, stim, true,  false);  % lungs only
[~,~,~,X_H] = build_fwd(P, P.maxsz, P.sig_heart, stim, false, true );  % heart only
Xh      = double(vh.meas(:)).';
d_lung  = X_L - Xh;
d_heart = X_H - Xh;
fprintf('  reference signatures: ||d_lung|| = %.4g, ||d_heart|| = %.4g, ', ...
        norm(d_lung), norm(d_heart));
fprintf('corr(d_lung,d_heart) = %.3f\n', pearson(d_lung,d_heart));

R1 = struct();
for m = 1:numel(D)
    C = D(m).comp;  K = size(C,1);
    nk = K - 1;                                  % last row is residual/mean
    rL = nan(1,nk);  rH = nan(1,nk);  sdv = nan(1,nk);
    for k = 1:nk
        rL(k) = abs(pearson(C(k,:), d_lung));
        rH(k) = abs(pearson(C(k,:), d_heart));
        w = ones(1,K); w(k) = 0;
        p = eval_GREIT_fig_merit(recon(imdl,vh,vi,C,w), tgt_heart);
        sdv(k) = p(4);
    end
    [~,k_ind] = max(rL ./ max(rH,eps));
    [~,k_rl]  = max(rL);
    [~,k_sd]  = min(sdv);
    fprintf('\n  %s  (%d components screened)\n', D(m).name, nk);
    fprintf('    %-24s %s\n','|corr| with lung signature', num2str(rL,'%8.3f'));
    fprintf('    %-24s %s\n','|corr| with heart signature',num2str(rH,'%8.3f'));
    fprintf('    %-24s %s\n','ratio lung/heart',           num2str(rL./max(rH,eps),'%8.3f'));
    fprintf('    %-24s %s\n','SD when removed',            num2str(sdv,'%8.3f'));
    fprintf('    max lung-corr #%d | max ratio #%d | min SD #%d  --> %s\n', ...
            k_rl, k_ind, k_sd, ternary(k_ind==k_sd,'AGREE','DISAGREE'));
    R1(m).name=D(m).name; R1(m).rL=rL; R1(m).rH=rH;
    R1(m).sd=sdv; R1(m).k_ind=k_ind; R1(m).k_rl=k_rl; R1(m).k_sd=k_sd;
end

%% =======================================================================
%  R6  the removed component alone
%  =======================================================================
fprintf('\n');
banner('R6  where the removed component puts its energy');
fprintf(['  Note: a component selected because its removal reveals the heart will\n' ...
        '  necessarily carry energy at the heart. Read against R1(a), which uses\n' ...
        '  measurement-domain ground truth and does not have that dependence.\n']);
fprintf('  %-8s %-6s %10s %10s %10s %10s\n', ...
        'method','comp','mean|.|lung','mean|.|heart','mean|.|bg','lung/heart');
R6 = struct();
for m = 1:numel(D)
    k  = D(m).key;
    ik = comp_image(imdl,vh,vi,D(m).comp,k);
    e  = ik.elem_data;
    el = mean(abs(e(roi_l)));  eh = mean(abs(e(roi_h)));  eb = mean(abs(e(roi_b)));
    fprintf('  %-8s #%-5d %10.4f %10.4f %10.4f %10.2f\n', ...
            D(m).name, k, el, eh, eb, el/eh);
    R6(m).name=D(m).name; R6(m).k=k; R6(m).lung=el; R6(m).heart=eh;
    R6(m).bg=eb; R6(m).img=e;
end

%% =======================================================================
%  R2  DAGAF window cap: full PE / SD / RNG triple
%  =======================================================================
fprintf('\n');
banner('R2  DAGAF window cap against the full triple');
caps = [2 3 4 6 8];
R2 = nan(numel(caps),4);   % cap, PE, SD, RNG
fprintf('  %5s %8s %10s %10s %10s\n','cap','window','PE','SD','RNG');
for c = 1:numel(caps)
    Dp = P.dagaf;  Dp.maskCap = caps(c);
    [cmp,~] = decomp_dagaf(X, Dp);
    K = size(cmp,1);
    w = ones(1,K); w(1) = 0;
    p = eval_GREIT_fig_merit(recon(imdl,vh,vi,cmp,w), tgt_heart);
    R2(c,:) = [caps(c) reshape(p(ridx),1,3)];
    fprintf('  %5d %8d %10.3f %10.3f %10.3f\n', ...
            caps(c), 2*caps(c)+1, p(2), p(4), p(5));
end

%% =======================================================================
%  R3  positive control: plain fixed linear high-pass filters
%      Two members of the class (Gaussian-subtraction, boxcar-subtraction),
%      two axes (concatenated channel index, injection index).
%  =======================================================================
fprintf('\n');
banner('R3  fixed linear high-pass controls');
ctrl = struct('kind',{},'axis',{},'L',{},'pe',{},'sd',{},'rng',{});
fprintf('  %-10s %-6s %4s %10s %10s %10s\n','filter','axis','L','PE','SD','RNG');
for kind = {'gauss','box'}
  for ax = {'chan','inj'}
    Ls = ternary(strcmp(ax{1},'chan'), [9 13], [5 7]);
    for L = Ls
        cmp = decomp_linhp(X, L, kind{1}, ax{1});
        err = norm(sum(cmp,1)-X)/norm(X);
        w = [0 1];                          % remove the high-pass row
        p = eval_GREIT_fig_merit(recon(imdl,vh,vi,cmp,w), tgt_heart);
        fprintf('  %-10s %-6s %4d %10.3f %10.3f %10.3f   (closure %.1e)\n', ...
                kind{1}, ax{1}, L, p(2), p(4), p(5), err);
        ctrl(end+1) = struct('kind',kind{1},'axis',ax{1},'L',L, ...
                             'pe',p(2),'sd',p(4),'rng',p(5)); %#ok<SAGROW>
    end
  end
end

%% =======================================================================
%  R4  variability across forward-mesh density
%      The GREIT reconstruction model is held fixed; only the forward mesh
%      changes. Component indices are pinned so that the sweep measures the
%      sensitivity of the metrics, not of the selection.
%  =======================================================================
fprintf('\n');
banner('R4  forward-mesh density sweep');
maxszs = [0.16 0.18 0.20 0.22 0.25];
methods = {'Original','PCA','EMD','DAGAF'};
R4 = nan(numel(maxszs), 4, 3);        % mesh x method x (PE SD RNG)
nel = nan(1,numel(maxszs));
for s = 1:numel(maxszs)
    [f2,vh2,vi2,X2] = build_fwd(P, maxszs(s), P.sig_heart, stim);
    nel(s) = size(f2.elems,1);
    D2 = decompose_all(X2, P);
    p = eval_GREIT_fig_merit(inv_solve(imdl,vh2,vi2), tgt_heart);
    R4(s,1,:) = reshape(p(ridx),1,1,3);
    for m = 1:numel(D2)
        C = D2(m).comp;  K = size(C,1);
        w = ones(1,K);  w(D2(m).key) = 0;
        if strcmp(D2(m).name,'PCA')
            keep = false(1,K); keep([P.pca_keep K]) = true;
            w = double(keep); w(D2(m).key) = 0;
        end
        p = eval_GREIT_fig_merit(recon(imdl,vh2,vi2,C,w), tgt_heart);
        R4(s,m+1,:) = reshape(p(ridx),1,1,3);
    end
    fprintf('  maxsz %.3f  (%6d elems)  done\n', maxszs(s), nel(s));
end
fprintf('\n  %-10s %-28s %-28s %-28s\n','method','PE  min .. max (spread)', ...
        'SD  min .. max (spread)','RNG min .. max (spread)');
for m = 1:4
    line = sprintf('  %-10s', methods{m});
    for q = 1:3
        v = R4(:,m,q);
        line = [line sprintf('%8.3f ..%8.3f (%6.3f)  ', min(v), max(v), max(v)-min(v))]; %#ok<AGROW>
    end
    fprintf('%s\n', line);
end

%% =======================================================================
%  R12  sigma_heart 1.5 against 2.0, full triple, against the published table
%  =======================================================================
fprintf('\n');
banner('R12  heart conductivity 1.5 vs 2.0 against the published triple');
pub = [0.04 0.98 38.77; 0.04 0.03 1.45; 0.04 0.02 0.77];   % Orig, PCA, EMD
sigs = [1.5 2.0];
R12 = nan(numel(sigs),4,3);
for s = 1:numel(sigs)
    [~,vh3,vi3,X3] = build_fwd(P, P.maxsz, sigs(s), stim);
    D3 = decompose_all(X3, P);
    p = eval_GREIT_fig_merit(inv_solve(imdl,vh3,vi3), tgt_heart);
    R12(s,1,:) = reshape(p(ridx),1,1,3);
    for m = 1:numel(D3)
        C = D3(m).comp;  K = size(C,1);
        if strcmp(D3(m).name,'PCA')
            keep = false(1,K); keep([P.pca_keep K]) = true;  w = double(keep);
        else
            w = ones(1,K);
        end
        w(D3(m).key) = 0;
        p = eval_GREIT_fig_merit(recon(imdl,vh3,vi3,C,w), tgt_heart);
        R12(s,m+1,:) = reshape(p(ridx),1,1,3);
    end
end
lbl = {'Original','PCA','EMD','DAGAF'};
for m = 1:4
    fprintf('  %-9s', lbl{m});
    for s = 1:2
        fprintf('  sigma %.1f: PE %6.3f SD %6.3f RNG %8.3f', ...
                sigs(s), R12(s,m,1), R12(s,m,2), R12(s,m,3));
    end
    if m <= 3
        fprintf('   | published PE %.2f SD %.2f RNG %.2f', pub(m,1),pub(m,2),pub(m,3));
    end
    fprintf('\n');
end
fprintf('\n  Distance to the published triple (|d| summed over SD and RNG,\n');
fprintf('  PE excluded because it does not reproduce at any mesh density):\n');
for m = 1:3
    d15 = abs(R12(1,m,2)-pub(m,2)) + abs(R12(1,m,3)-pub(m,3));
    d20 = abs(R12(2,m,2)-pub(m,2)) + abs(R12(2,m,3)-pub(m,3));
    fprintf('    %-9s sigma1.5 %8.3f   sigma2.0 %8.3f   -> %s\n', ...
            lbl{m}, d15, d20, ternary(d15<d20,'1.5 closer','2.0 closer'));
end

%% ---- save --------------------------------------------------------------
save('revision_r1_results.mat','R1','R2','R4','R6','R12','ctrl', ...
     'caps','maxszs','nel','sigs','pub','P');
fprintf('\nTotal elapsed %.1f min\n', toc(t_all)/60);
fprintf('wrote revision_r1_results.txt and revision_r1_results.mat\n');
diary off

%% ========================================================================
function banner(s)
fprintf('=====================================================================\n');
fprintf(' %s\n', s);
fprintf('=====================================================================\n');
end
% -------------------------------------------------------------------------
function out = ternary(c,a,b), if c, out=a; else, out=b; end, end
% -------------------------------------------------------------------------
function r = pearson(x,y)
x = x(:) - mean(x(:));  y = y(:) - mean(y(:));
r = (x.'*y) / (norm(x)*norm(y));
end
% -------------------------------------------------------------------------
function [fmdl,vh,vi,X] = build_fwd(P, maxsz, sig_heart, stim, useLung, useHeart)
% useLung/useHeart default to true: the full phantom. Setting one false gives
% the isolated forward solution used for the ground-truth signatures of R1(a).
if nargin < 5, useLung  = true; end
if nargin < 6, useHeart = true; end
lungs = sprintf('solid obj1 = sphere(1,0.1,1;%.2f) or sphere(-1,0.1,1;%.2f);', ...
                P.r_lung, P.r_lung);
heart = sprintf('solid obj2 = sphere(0,0.1,1;%.2f);', P.r_heart);
[fmdl, midx] = ng_mk_cyl_models([2,2,maxsz],[16,1],[0.1], ...
                                {'obj1','obj2',[lungs heart]});
fmdl.stimulation = stim;
img = mk_image(fmdl,1);
vh  = fwd_solve(img);
if useLung,  img.elem_data(midx{2}) = P.sig_lung;  end
if useHeart, img.elem_data(midx{3}) = sig_heart;   end
vi  = fwd_solve(img);
X   = double(vi.meas(:)).';
end
% -------------------------------------------------------------------------
function D = decompose_all(X, P)
D = struct('name',{},'comp',{},'key',{});
[cp,~] = decomp_pca(X);
D(end+1) = struct('name','PCA','comp',cp,'key',P.pca_lung_key);
[ce,~] = decomp_emd(X, P.emd_nIMF);
D(end+1) = struct('name','EMD','comp',ce,'key',P.emd_lung_key);
[cd_,~] = decomp_dagaf(X, P.dagaf);
D(end+1) = struct('name','DAGAF','comp',cd_,'key',1);
end
% -------------------------------------------------------------------------
function [comp, note] = decomp_pca(X)
A  = reshape(X(:),13,16);
mu = mean(A,1);
Ac = A - repmat(mu,13,1);
[V,Dg] = eig(cov(Ac));
[~,ord] = sort(diag(Dg),'descend');
V  = V(:,ord);
PC = Ac*V;
K  = size(V,2);
comp = zeros(K+1,numel(X));
for k = 1:K, comp(k,:) = reshape(PC(:,k)*V(:,k)',1,[]); end
comp(K+1,:) = reshape(repmat(mu,13,1),1,[]);
note = sprintf('%d components + mean row',K);
end
% -------------------------------------------------------------------------
function [comp, note] = decomp_emd(X, nIMF)
[imf,res] = emd_rcada(X);
comp = [imf; res];
note = sprintf('%d IMFs + residual',size(imf,1));
end
% -------------------------------------------------------------------------
function [comp, note] = decomp_dagaf(X, Dp)
[imf,res,info] = dagaf2(X, Dp.nIMF, Dp.chi, Dp.eType, Dp.eType, ...
                        Dp.TH1, Dp.TH2, Dp.alpha, Dp.maskCap);
comp = [imf; res];
note = sprintf('%d IMFs + residual, M = [%s]', size(imf,1), num2str(info.M));
end
% -------------------------------------------------------------------------
function comp = decomp_linhp(X, L, kind, axis_)
% Two-row decomposition [high-pass; low-pass] summing exactly to X, where the
% low-pass is a plain fixed linear convolution with a normalised window and
% periodic extension. This is the positive control for the claim that the
% ringing penalty follows from removing a fixed linear high-pass band.
X = reshape(double(X),1,[]);
switch axis_
    case 'chan'
        lo = linsmooth(X, L, kind);
    case 'inj'
        A = reshape(X,13,16);
        lo = zeros(13,16);
        for i = 1:13, lo(i,:) = linsmooth(A(i,:), L, kind); end
        lo = reshape(lo,1,208);
end
comp = [X - lo; lo];
end
% -------------------------------------------------------------------------
function y = linsmooth(x, L, kind)
x = x(:).';  n = numel(x);
if mod(L,2)==0, L = L+1; end
h = (L-1)/2;
switch kind
    case 'gauss'
        alpha = 3;  k = -h:h;
        w = exp(-0.5*(alpha*k/h).^2);
    case 'box'
        w = ones(1,L);
end
w  = w/sum(w);
xp = [x(n-h+1:n), x, x(1:h)];        % periodic extension, as DAGAF eType 'p'
y  = conv(xp, w, 'valid');
end
% -------------------------------------------------------------------------
function img = recon(imdl, vh, vi, comp, w)
v = vi;  v.meas = (w(:).'*comp).';
img = inv_solve(imdl, vh, v);
img.calc_colours.ref_level = 0;
end
% -------------------------------------------------------------------------
function imgk = comp_image(imdl, vh, vi, comp, k)
% The contribution of component k to the reconstruction. One-step difference
% reconstruction is linear in the measurement difference, so this is exact.
K  = size(comp,1);
w0 = ones(1,K);  wk = w0;  wk(k) = 0;
a  = recon(imdl,vh,vi,comp,w0);
b  = recon(imdl,vh,vi,comp,wk);
imgk = a;  imgk.elem_data = a.elem_data - b.elem_data;
end
