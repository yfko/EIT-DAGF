%% EIT_REVISION_R4_LUNG
% Completes roadmap item R4 for the lung target: figure-of-merit ranges
% across forward-mesh density, with the negative-contrast correction of
% Section 2.6 applied and the reconstruction model held fixed.
%
% Writes revision_r4_lung.txt and revision_r4_lung.mat.

clc; clear; close all;
diary off;
if exist('revision_r4_lung.txt','file'), delete('revision_r4_lung.txt'); end
diary revision_r4_lung.txt
t_all = tic;

P.n_elecs   = 16;
P.r_heart   = 0.3;
P.sig_lung  = 0.5;
P.sig_heart = 2.0;
P.w_boost   = 1.001;
P.pca_keep     = 1:12;
P.pca_lung_key = 4;
P.emd_lung_key = 1;
P.dagaf = struct('nIMF',6,'chi',2,'alpha',3,'eType','p', ...
                 'maskCap',6,'TH1',40,'TH2',1e-8);

here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(here,'rcada_emd'));
addpath(fullfile(here,'github_release','src'));
stim = mk_stim_patterns(P.n_elecs,1,[0,1],[0,1],{'no_meas_current'},1);
ridx = [2 4 5];

fprintf('Building GREIT model ...\n');
opt   = struct('imgsz',[64 64],'distr',0,'Nsim',1000, ...
               'target_size',0.02,'noise_figure',0.5);
fmdlG = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1]);
fmdlG.stimulation = stim;
fmdlG = mdl_normalize(fmdlG,1);
imdl  = mk_GREIT_model(mk_image(fmdlG,1),0.2,[],opt);
ctrT  = interp_mesh(imdl.rec_model);
ctrR  = (ctrT(1:2:end,:) + ctrT(2:2:end,:))/2;
mask_right = ctrR(:,1) > 0;          % keep one lung, in image coordinates

maxszs  = [0.16 0.18 0.20 0.22 0.25];
radii   = [0.6 0.3];
methods = {'Original','PCA','EMD','DAGAF'};
R = nan(numel(radii), numel(maxszs), 4, 3);

for ri = 1:numel(radii)
  P.r_lung = radii(ri);
  tgt_lung = [1; 0.1; 1; P.r_lung];
  lungfom  = @(im) eval_GREIT_fig_merit(negmask(im, mask_right), tgt_lung);
  for s = 1:numel(maxszs)
    [~,vh,vi,X] = build_fwd(P, maxszs(s), stim);
    D = decompose_all(X, P);
    p = lungfom(inv_solve_ref(imdl,vh,vi));
    R(ri,s,1,:) = reshape(p(ridx),1,1,1,3);
    for m = 1:numel(D)
        C = D(m).comp;  K = size(C,1);
        if strcmp(D(m).name,'PCA')
            keep = false(1,K); keep([P.pca_keep K]) = true; w = double(keep);
        else
            w = ones(1,K);
        end
        w(D(m).key) = P.w_boost;                  % lung imaging: boost it
        p = lungfom(recon(imdl,vh,vi,C,w));
        R(ri,s,m+1,:) = reshape(p(ridx),1,1,1,3);
    end
    fprintf('  r_lung %.1f  maxsz %.3f  done\n', P.r_lung, maxszs(s));
  end
end

names = {'PE','SD','RNG'};
for ri = 1:numel(radii)
  fprintf('\n===== r_lung = %.1f =====\n', radii(ri));
  for q = 1:3
    fprintf('  %-4s', names{q});
    for m = 1:4
        v = squeeze(R(ri,:,m,q));
        fprintf('  %-9s %6.3f (%6.3f-%6.3f)', methods{m}, v(3), min(v), max(v));
    end
    fprintf('\n');
  end
end

save('revision_r4_lung.mat','R','maxszs','radii','methods');
fprintf('\nElapsed %.1f min\n', toc(t_all)/60);
diary off

%% ========================================================================
function [fmdl,vh,vi,X] = build_fwd(P, maxsz, stim)
lungs = sprintf('solid obj1 = sphere(1,0.1,1;%.2f) or sphere(-1,0.1,1;%.2f);', ...
                P.r_lung, P.r_lung);
heart = sprintf('solid obj2 = sphere(0,0.1,1;%.2f);', P.r_heart);
[fmdl, midx] = ng_mk_cyl_models([2,2,maxsz],[16,1],[0.1], ...
                                {'obj1','obj2',[lungs heart]});
fmdl.stimulation = stim;
img = mk_image(fmdl,1);
vh  = fwd_solve(img);
img.elem_data(midx{2}) = P.sig_lung;
img.elem_data(midx{3}) = P.sig_heart;
vi  = fwd_solve(img);
X   = double(vi.meas(:)).';
end
% -------------------------------------------------------------------------
function D = decompose_all(X, P)
D = struct('name',{},'comp',{},'key',{});
A  = reshape(X(:),13,16); mu = mean(A,1); Ac = A - repmat(mu,13,1);
[V,Dg] = eig(cov(Ac)); [~,ord] = sort(diag(Dg),'descend'); V = V(:,ord);
PC = Ac*V; K = size(V,2);
cp = zeros(K+1,numel(X));
for k = 1:K, cp(k,:) = reshape(PC(:,k)*V(:,k)',1,[]); end
cp(K+1,:) = reshape(repmat(mu,13,1),1,[]);
D(end+1) = struct('name','PCA','comp',cp,'key',P.pca_lung_key);
[imf,res] = emd_rcada(X);
D(end+1) = struct('name','EMD','comp',[imf; res],'key',P.emd_lung_key);
Dp = P.dagaf;
[imf2,res2] = dagaf2(X, Dp.nIMF, Dp.chi, Dp.eType, Dp.eType, Dp.TH1, Dp.TH2, Dp.alpha, Dp.maskCap);
D(end+1) = struct('name','DAGAF','comp',[imf2; res2],'key',1);
end
% -------------------------------------------------------------------------
function img = inv_solve_ref(imdl, vh, vi)
img = inv_solve(imdl, vh, vi); img.calc_colours.ref_level = 0;
end
% -------------------------------------------------------------------------
function img = recon(imdl, vh, vi, comp, w)
v = vi; v.meas = (w(:).'*comp).';
img = inv_solve(imdl, vh, v); img.calc_colours.ref_level = 0;
end
% -------------------------------------------------------------------------
function im = negmask(im, keep)
im.elem_data = -im.elem_data;
im.elem_data(~keep) = 0;
end
