%% EIT_COMPARE_TABLE1
% Head-to-head comparison of three pre-processing methods applied to the raw
% 208 EIT voltages before reconstruction, reported in the format of Table 1 of
%
%   K.-S. Cheng, P.-L. Su, Y.-F. Ko, "Separation of Heart and Lung-related
%   Signals in Electrical Impedance Tomography Using Empirical Mode
%   Decomposition", Curr Med Imaging 2022;18(13):1396-1415.
%
%   PCA    - as in that paper (208 -> 13x16, 16 components, PC4 = lung key)
%   EMD    - as in that paper (6 IMFs, IMF1 = lung key)
%   DAGAF  - Lin, Tan & Tian, Biomed Signal Process Control 71:103104, 2022
%
% Geometry, conductivities and GREIT settings are taken verbatim from
% "Error quantification_GREIT Test Parameters.m".
%
% UNIFYING IDEA
%   Every method is reduced to the same interface: a K x 208 component matrix
%   whose rows sum exactly to the measured vector. Heart imaging = drop the
%   lung-dominant row; lung imaging = scale it up slightly. Because all three
%   then differ only in *which* decomposition produced the rows, the
%   comparison is apples-to-apples.
%
% Requires: EIDORS >= 3.8 on the path, dagaf2.m in this folder.
%           EMD is optional (MATLAB Signal Processing Toolbox emd(), or
%           Rilling's emd.m). If absent, the EMD column reports NaN.
%
% Outputs:  printed Table 1, plus table1_comparison.csv and a figure.

clc; clear; close all;
% startup;   % <- run EIDORS startup first if not already on the path

%% ---- parameters --------------------------------------------------------
P.n_elecs   = 16;
P.r_lung    = 0.6;      % Table 1 in the paper is for r_lung = 0.6
P.r_heart   = 0.3;
P.sig_lung  = 0.5;
P.sig_heart = 2.0;    % value stated in the 2022 Methods. Table 1 there was in fact
                       % produced with 1.5 (see the sigma check in the reproducibility section).
P.cut_row   = 31;
P.xaxis     = -2:4/63:2;

P.w_boost   = 1.001;    % lung imaging: scale the lung-key component
                        % (original script used 1.001; paper text says 1.01)

% --- PCA (as described in the CMI paper) ---
P.pca_keep     = 1:12;  % Fig.5: components 1-12 are needed for a clean image
P.pca_lung_key = 4;     % PC4

% --- EMD ---
P.emd_nIMF     = 6;
P.emd_lung_key = 1;     % IMF1

% --- DAGAF ---
P.dagaf = struct('nIMF',6, 'chi',2, 'alpha',3, 'eType','p', ...
                 'maskCap',6, 'TH1',40, 'TH2',1e-8);   % 6 = one injection block
P.dagaf_lung_key = [];  % [] = decide by screening; or pin it to an integer

P.autoscreen = true;    % also screen every method and report what it picks

addpath(fullfile(fileparts(mfilename('fullpath')),'rcada_emd'));   % RCADA eemd
stim = mk_stim_patterns(P.n_elecs, 1, [0,1], [0,1], {'no_meas_current'}, 1);

%% ---- 1. phantom --------------------------------------------------------
fprintf('Building phantom (r_lung = %.1f)...\n', P.r_lung);
lungs = sprintf('solid obj1 = sphere(1,0.1,1;%.2f) or sphere(-1,0.1,1;%.2f);', ...
                P.r_lung, P.r_lung);
heart = sprintf('solid obj2 = sphere(0,0.1,1;%.2f);', P.r_heart);
[fmdl, midx] = ng_mk_cyl_models([2,2,0.2], [16,1], [0.1], {'obj1','obj2',[lungs heart]});
fmdl.stimulation = stim;

img = mk_image(fmdl, 1);
vh  = fwd_solve(img);
img.elem_data(midx{2}) = P.sig_lung;
img.elem_data(midx{3}) = P.sig_heart;
vi  = fwd_solve(img);
sim_slice = calc_slices(img);
sim_cut   = sim_slice(P.cut_row,:);

X = double(vi.meas(:)).';
N = numel(X);
fprintf('  %d measurements per frame\n', N);

%% ---- 2. GREIT reconstruction model ------------------------------------
fprintf('Building GREIT model...\n');
opt   = struct('imgsz',[64 64], 'distr',0, 'Nsim',1000, ...
               'target_size',0.02, 'noise_figure',0.5);
fmdlG = ng_mk_cyl_models([2,2,0.2], [16,1], [0.1]);
fmdlG.stimulation = stim;
fmdlG = mdl_normalize(fmdlG, 1);
imdl  = mk_GREIT_model(mk_image(fmdlG,1), 0.2, [], opt);

% GREIT's rec_model splits each square pixel into TWO triangles, but
% img.elem_data carries one value per PIXEL (here 6456 triangles -> 3228
% pixels). Averaging the paired triangle centres gives coordinates that line
% up with elem_data. Verified: paired centres are 0.0295 apart vs a 0.0625
% pixel pitch, so the 1:2:end / 2:2:end pairing is correct.
ctrT  = interp_mesh(imdl.rec_model);
ctrR  = (ctrT(1:2:end,:) + ctrT(2:2:end,:)) / 2;

%% ---- 3. decompose with each method ------------------------------------
M = struct('name',{},'comp',{},'key',{},'keep',{},'note',{});

% --- PCA ---
[cp, note] = decomp_pca(X);
M(end+1) = struct('name','PCA', 'comp',cp, 'key',P.pca_lung_key, ...
                  'keep',[P.pca_keep, size(cp,1)], 'note',note);   % last row = mean

% --- EMD ---
[ce, note] = decomp_emd(X, P.emd_nIMF);
if isempty(ce)
    M(end+1) = struct('name','EMD','comp',[],'key',[],'keep',[],'note',note);
else
    M(end+1) = struct('name','EMD', 'comp',ce, 'key',P.emd_lung_key, ...
                      'keep',1:size(ce,1), 'note',note);
end

% --- DAGAF ---
[cd_, note] = decomp_dagaf(X, P.dagaf);
M(end+1) = struct('name','DAGAF', 'comp',cd_, 'key',P.dagaf_lung_key, ...
                  'keep',1:size(cd_,1), 'note',note);

fprintf('\n--- decompositions ---\n');
for m = 1:numel(M)
    fprintf('  %-6s %s\n', M(m).name, M(m).note);
    if ~isempty(M(m).comp)
        err = norm(sum(M(m).comp,1) - X) / norm(X);
        fprintf('         completeness ||sum(comp)-X||/||X|| = %.2e\n', err);
    end
end

%% ---- 4. ROIs and baseline ----------------------------------------------
roi_h = sum((ctrR(:,1:2) - [0  0.1]).^2, 2) < P.r_heart^2;
roi_l = sum((ctrR(:,1:2) - [1  0.1]).^2, 2) < P.r_lung^2 | ...
        sum((ctrR(:,1:2) - [-1 0.1]).^2, 2) < P.r_lung^2;
% Keep only the x > 0 half, so eval_GREIT_fig_merit sees ONE lung.
% NOTE: the original script did this with elem_data(3228/2:end) = 0. Pixels
% are ordered by y, so that masks the TOP half (y > 0), leaving the lower
% halves of BOTH lungs rather than isolating one -- verified empirically.
mask_right = ctrR(:,1) > 0;

img_orig = inv_solve(imdl, vh, vi);
img_orig.calc_colours.ref_level = 0;

%% ---- 5. per-method screening -------------------------------------------
tgt_heart_scr = [0; 0.1; 1; P.r_heart];
if P.autoscreen
    fprintf('\n--- component screening by GREIT shape deformation ---\n');
    for m = 1:numel(M)
        if isempty(M(m).comp), continue; end
        K = size(M(m).comp,1);
        % Select by the actual GREIT figure of merit: the lung-dominant
        % component is the one whose removal minimises shape deformation of
        % the heart target. (An earlier ROI-ratio heuristic was discarded --
        % it rewarded destroying the image, so it picked PC1 over PC4.)
        nk = K - 1;                      % last row is the residual/mean
        sd = inf(1,nk);
        for k = 1:nk
            w = ones(1,K); w(k) = 0;
            pk = eval_GREIT_fig_merit(recon(imdl,vh,vi,M(m).comp,w), tgt_heart_scr);
            sd(k) = pk(4);
        end
        [~, kbest] = min(sd);
        fprintf('       SD by dropped component: %s\n', num2str(sd,'%.3f  '));
        fprintf('  %-6s screened lung-key = #%d', M(m).name, kbest);
        if ~isempty(M(m).key)
            fprintf('   (paper/config uses #%d)%s\n', M(m).key, ...
                    repmat(' <-- MISMATCH', 1, kbest ~= M(m).key));
        else
            fprintf('   (unset -> using screened value)\n');
            M(m).key = kbest;
        end
        M(m).screened = kbest;   %#ok<SAGROW>
    end
end

%% ---- 6. reconstruct heart and lung images for every method ------------
imgs = {img_orig};  names = {'Original'};
for m = 1:numel(M)
    if isempty(M(m).comp)
        imgs{end+1} = [];  names{end+1} = [M(m).name '-based'];  %#ok<SAGROW>
        continue
    end
    K = size(M(m).comp,1);
    if isempty(M(m).key)
        warning('%s: no lung-key set and autoscreen is off; defaulting to #1', M(m).name);
        M(m).key = 1;
    end
    keep = false(1,K); keep(M(m).keep) = true;

    w_heart = double(keep);  w_heart(M(m).key) = 0;              % drop lung key
    w_lung  = double(keep);  w_lung(M(m).key)  = P.w_boost;      % boost it

    imgs{end+1}  = recon(imdl,vh,vi,M(m).comp,w_heart);          %#ok<SAGROW>
    names{end+1} = [M(m).name '-based'];                         %#ok<SAGROW>
    IL(m).img    = recon(imdl,vh,vi,M(m).comp,w_lung);           %#ok<SAGROW>
end

%% ---- 7. TABLE 1 --------------------------------------------------------
tgt_heart = [0; 0.1; 1; P.r_heart];
tgt_lung  = [1; 0.1; 1; P.r_lung];

rows = {'PE','SD','RNG'};   ridx = [2 4 5];      % eval_GREIT: AR PE RES SD RNG
T = nan(3, numel(imgs));
for j = 1:numel(imgs)
    if isempty(imgs{j}), continue; end
    p = eval_GREIT_fig_merit(imgs{j}, tgt_heart);
    T(:,j) = p(ridx);
end

fprintf('\n');
fprintf('=========================================================================\n');
fprintf(' Table 1.  PE, SD and RNG for HEART imaging (r_lung = %.1f, r_heart = %.1f)\n', ...
        P.r_lung, P.r_heart);
fprintf('=========================================================================\n');
fprintf('%-8s', 'Heart');
fprintf('%14s', names{:});  fprintf('\n');
for i = 1:3
    fprintf('%-8s', rows{i});
    for j = 1:numel(imgs), fprintf('%14.3f', T(i,j)); end
    fprintf('\n');
end
fprintf('-------------------------------------------------------------------------\n');
fprintf(' CMI 2022 Table 1 (published):  Original  PE 0.04  SD 0.98  RNG 38.77\n');
fprintf('                                PCA-based PE 0.04  SD 0.03  RNG  1.45\n');
fprintf('                                EMD-based PE 0.04  SD 0.02  RNG  0.77\n');
fprintf('=========================================================================\n');

% --- same table for LUNG imaging ---
TL = nan(3, numel(M)+1);
% The lung is a NEGATIVE-contrast target (sigma 0.5 < background 1), but
% eval_GREIT_fig_merit locates a positive blob. Without negating the image
% first it reports SD = 1.000 for every method -- which is what the earlier
% version of this table did. Negate, then keep one lung.
lungfom = @(im) eval_GREIT_fig_merit(negmask(im, mask_right), tgt_lung);
p = lungfom(img_orig);  TL(:,1) = p(ridx);
for m = 1:numel(M)
    if isempty(M(m).comp), continue; end
    p = lungfom(IL(m).img);  TL(:,m+1) = p(ridx);
end
fprintf('\n Table 1b. PE, SD and RNG for LUNG imaging (one lung masked out)\n');
fprintf('%-8s', 'Lung');  fprintf('%14s', names{:});  fprintf('\n');
for i = 1:3
    fprintf('%-8s', rows{i});
    for j = 1:size(TL,2), fprintf('%14.3f', TL(i,j)); end
    fprintf('\n');
end
fprintf('=========================================================================\n');

%% ---- 8. CSV ------------------------------------------------------------
fid = fopen('table1_comparison.csv','w');
fprintf(fid, 'target,metric,%s\n', strjoin(names, ','));
for i = 1:3
    fprintf(fid, 'heart,%s', rows{i});
    fprintf(fid, ',%.4f', T(i,:));   fprintf(fid, '\n');
end
for i = 1:3
    fprintf(fid, 'lung,%s', rows{i});
    fprintf(fid, ',%.4f', TL(i,:));  fprintf(fid, '\n');
end
fclose(fid);
fprintf('\nwrote table1_comparison.csv\n');

%% ---- 9. figure ---------------------------------------------------------
n = numel(imgs);
fprintf('\n--- centre-cut profiles on the simulated scale (phantom: 0.5 .. 2.0) ---\n');
figure('Color','w','Name','PCA / EMD / DAGAF comparison','Position',[80 80 320*n 520]);
for j = 1:n
    subplot(2,n,j);
    if isempty(imgs{j}), axis off; title([names{j} ' (n/a)']); continue; end
    show_slices(imgs{j}); title(names{j});
    subplot(2,n,j+n);
    sl = calc_slices(imgs{j});
    prof = normprof(sl(P.cut_row,:), sim_cut);
    fprintf('  %-12s centre-cut peak %.2f, flanks %.2f / %.2f\n', names{j}, ...
            max(prof,[],'omitnan'), prof(16), prof(49));
    % same styling as the three-state figure: phantom in solid black first,
    % reconstruction in blue over it, one fixed axis in every panel
    plot(P.xaxis, sim_cut, 'k-', 'LineWidth',1.3); hold on
    plot(P.xaxis, prof, 'b-', 'LineWidth',1.7);
    grid on; axis square; xlabel('x'); ylabel('\sigma');
    ylim([0 2.1]); yticks(0:0.5:2);
    if j == 1
        legend('simulated','reconstructed','Location','south'); legend boxoff
    end
end

exportgraphics(gcf,'fig_compare_profiles.png','Resolution',150);
fprintf('wrote fig_compare_profiles.png\n');

%% ========================================================================
function [comp, note] = decomp_pca(X)
% CMI paper 2.2: reshape 208 -> 13 x 16 (column = injection pair),
% PCA over the 16 injections -> 16 components. Rows of comp are the rank-1
% reconstructions; the last row carries the mean so the rows sum to X.
A  = reshape(X(:), 13, 16);
mu = mean(A, 1);
Ac = A - repmat(mu, 13, 1);
Cv = cov(Ac);                       % 16 x 16
[V, D] = eig(Cv);
[~, ord] = sort(diag(D), 'descend');
V  = V(:, ord);
PC = Ac * V;                        % 13 x 16 scores
K  = size(V,2);
comp = zeros(K+1, numel(X));
for k = 1:K
    comp(k,:) = reshape(PC(:,k) * V(:,k)', 1, []);
end
comp(K+1,:) = reshape(repmat(mu,13,1), 1, []);
note = sprintf('%d components + mean row', K);
end
% -------------------------------------------------------------------------
function [comp, note] = decomp_emd(X, nIMF)
% Works with both MATLAB's built-in emd() (returns N x K) and Rilling's
% emd.m (returns K x N). The residual is appended so the rows sum to X.
comp = [];  note = '';
if exist('eemd','file') == 2
    % The authentic RCADA/NCU package (recovered from the Internet Archive,
    % see emd_rcada.m). This is the code the CMI 2022 paper used.
    [imf, res] = emd_rcada(X);
    comp = [imf; res];
    note = sprintf('%d IMFs + residual (RCADA eemd, Nstd=0 NE=1)', size(imf,1));
    return
end
if exist('emd','file') == 0
    [imf, res] = emd_basic(X, nIMF);
    comp = [imf; res];
    note = sprintf('%d IMFs + residual (emd_basic fallback)', size(imf,1));
    return
end
try
    r = emd(X(:), 'MaxNumIMF', nIMF);
catch
    try
        r = emd(X(:));
    catch ME
        note = ['emd() failed: ' ME.message];  return
    end
end
if size(r,2) == numel(X), imf = r; else, imf = r.'; end
if size(imf,1) > nIMF                       % fold the extras into the last row
    imf(nIMF,:) = sum(imf(nIMF:end,:), 1);
    imf(nIMF+1:end,:) = [];
end
res  = X(:).' - sum(imf,1);
comp = [imf; res];
note = sprintf('%d IMFs + residual row', size(imf,1));
end
% -------------------------------------------------------------------------
function [comp, note] = decomp_dagaf(X, D)
[imf, res, info] = dagaf2(X, D.nIMF, D.chi, D.eType, D.eType, ...
                          D.TH1, D.TH2, D.alpha, D.maskCap);
comp = [imf; res];
note = sprintf('%d IMFs + residual row, M = [%s], Ne = [%s] (%s)', ...
               size(imf,1), num2str(info.M), num2str(info.Ne), info.stop);
end
% -------------------------------------------------------------------------
function img = recon(imdl, vh, vi, comp, w)
v = vi;  v.meas = (w(:).' * comp).';
img = inv_solve(imdl, vh, v);
img.calc_colours.ref_level = 0;
end
% -------------------------------------------------------------------------
function s = roi_vals(img, roi_h, roi_l)
s.h = mean(img.elem_data(roi_h));
s.l = mean(img.elem_data(roi_l));
end
% -------------------------------------------------------------------------
function y = normprof(x, sim)
% Put the reconstructed profile on the simulated conductivity scale so the two
% can be overlaid: remove DC, scale to the simulated range, re-centre on it.
% The 2022 scripts divided by the reconstruction's own range and added 1, which
% only lines up when sigma_heart = 1.5 makes the simulated range equal to 1.
x = x - mean(x(~isnan(x)));
r = max(x(:),[],'omitnan') - min(x(:),[],'omitnan');
sr = max(sim(:),[],'omitnan') - min(sim(:),[],'omitnan');
y = x/r*sr + mean(sim(~isnan(sim)));
end
% -------------------------------------------------------------------------
function im = negmask(im, keep)
% Flip sign (negative-contrast target) and zero everything outside `keep`.
im.elem_data = -im.elem_data;
im.elem_data(~keep) = 0;
end
