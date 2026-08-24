%% EIT_DAGAF_IMAGES  Render the DAGAF-processed EIT reconstructions.
% Produces:
%   fig_dagaf_separation.png  - phantom / original / PCA heart / DAGAF heart / DAGAF lung
%   fig_dagaf_components.png  - the DAGAF decomposition of the 208-voltage vector
% Run setup_eit_env first.

clc; close all;
P.r_lung = 0.6; P.r_heart = 0.3; P.sig_lung = 0.5; P.sig_heart = 2.0;
P.k_lung_dagaf = 1;    % established by the FoM sweep
P.k_lung_pca   = 4;    % as in Cheng, Su & Ko 2022
P.w_boost      = 1.05;   % 1.001 (original script) changes the image by <1%;
                         % 1.05 is where the heart stops being distinguishable
xax = -2:4/63:2;

stim = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);

%% phantom
lungs = sprintf('solid obj1 = sphere(1,0.1,1;%.2f) or sphere(-1,0.1,1;%.2f);',P.r_lung,P.r_lung);
heart = sprintf('solid obj2 = sphere(0,0.1,1;%.2f);',P.r_heart);
[fmdl,midx] = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1],{'obj1','obj2',[lungs heart]});
fmdl.stimulation = stim;
imgp = mk_image(fmdl,1); vh = fwd_solve(imgp);
imgp.elem_data(midx{2}) = P.sig_lung; imgp.elem_data(midx{3}) = P.sig_heart;
vi = fwd_solve(imgp);
truth = calc_slices(imgp);
X = double(vi.meas(:)).';

%% GREIT
opt = struct('imgsz',[64 64],'distr',0,'Nsim',1000,'target_size',0.02,'noise_figure',0.5);
fg = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1]); fg.stimulation = stim; fg = mdl_normalize(fg,1);
imdl = mk_GREIT_model(mk_image(fg,1),0.2,[],opt);

%% decompose
A = reshape(X(:),13,16); mu = mean(A,1); Ac = A-repmat(mu,13,1);
[V,L] = eig(cov(Ac)); [~,o] = sort(diag(L),'descend'); V = V(:,o); PCs = Ac*V;
Cp = zeros(17,208);
for k=1:16, Cp(k,:) = reshape(PCs(:,k)*V(:,k)',1,[]); end
Cp(17,:) = reshape(repmat(mu,13,1),1,[]);
% nIMF=2: verified that PE/SD/RNG are bit-identical for nIMF = 2,3,4,6,8,
% because dropping IMF1 leaves X - IMF1 however the rest is subdivided.
[di,dr] = dagaf2(X,2,2,'p','p',40,1e-8,3,6);
Cd = [di;dr];

[ei,er] = emd_rcada(X);          % RCADA eemd, same input as the others
Ce = [ei;er];

rc = @(C,w) inv_solve(imdl, vh, setfield(vi,'meas',(w*C).'));   %#ok<SFLD>
wP = ones(1,17);         wP(P.k_lung_pca)   = 0;
wE = ones(1,size(Ce,1)); wE(1)              = 0;
wD = ones(1,size(Cd,1)); wD(P.k_lung_dagaf) = 0;

% Same four columns as Figure 2, so the two figures can be read against each
% other: the unprocessed frame and the three heart reconstructions.
R = {calc_slices(inv_solve(imdl,vh,vi)), calc_slices(rc(Cp,wP)), ...
     calc_slices(rc(Ce,wE)), calc_slices(rc(Cd,wD))};
T = {'Original (mixed)','PCA-based heart (drop PC4)', ...
     'EMD-based heart (drop IMF1)','DAGAF-based heart (drop IMF1)'};

%% ---- figure 1: separation -------------------------------------------------
cm = bwr(256);
f = figure('Color','w','Position',[0 0 1500 640],'Visible','off');
% Panel 1: the Netgen mesh itself. The triangles, nodes and region indices
% come straight from fmdl (ng_mk_cyl_models -> Netgen); only the renderer is
% ours, because EIDORS' show_fem draws unclipped slivers across the axes --
% verified separately that all its patches are in-bounds, so the streaks are
% a plotting artefact and not the geometry.
tmpd = tempname; mkdir(tmpd); pmesh = fullfile(tmpd,'mesh.png');
hm = figure('Color','w','Position',[0 0 760 760],'Visible','off');
fmv = imgp.fwd_model;  Nv = fmv.nodes;  sg = imgp.elem_data;
% fmv.boundary is the boundary EIDORS stored when Netgen built the mesh.
% Recomputing it with find_boundary() returns fewer faces (5394 vs 5988) and
% draws long stray edges across the axes -- use the stored one.
trisurf(fmv.boundary, Nv(:,1),Nv(:,2),Nv(:,3), ...
        'FaceColor',[0.88 0.89 0.92],'FaceAlpha',0.30, ...
        'EdgeColor',[0.35 0.38 0.44],'LineWidth',0.25); hold on
trisurf(find_boundary(fmv.elems(abs(sg-0.5)<1e-9,:)), Nv(:,1),Nv(:,2),Nv(:,3), ...
        'FaceColor',[0.20 0.32 0.80],'FaceAlpha',0.92,'EdgeColor','none');
trisurf(find_boundary(fmv.elems(abs(sg-P.sig_heart)<1e-9,:)), Nv(:,1),Nv(:,2),Nv(:,3), ...
        'FaceColor',[0.85 0.18 0.18],'FaceAlpha',0.95,'EdgeColor','none');
for e = 1:numel(fmv.electrode)
    ce = mean(Nv(fmv.electrode(e).nodes,:),1);
    plot3(ce(1),ce(2),ce(3),'o','MarkerSize',7, ...
          'MarkerFaceColor',[0.10 0.62 0.20],'MarkerEdgeColor','k');
end
view(2); axis equal off; xlim([-2.15 2.15]); ylim([-2.15 2.15]);
exportgraphics(hm,pmesh,'Resolution',220,'ContentType','image'); close(hm);

subplot(2,5,1); imshow(imread(pmesh));
title({'FEM mesh (Netgen)',sprintf('\\sigma: heart %.1f, lung %.1f',P.sig_heart,P.sig_lung)},'FontSize',10);
subplot(2,5,6); plot(xax,truth(31,:),'k','LineWidth',1.5); axis square; grid on;
xlabel('x'); ylabel('\sigma'); ylim([0 2.2]); title('phantom centre cut','FontSize',9);
% Absolute Delta-sigma is not comparable across methods: each modification
% changes the overall magnitude of the measured difference, so a shared axis
% would hide the weaker panels while independent axes misrepresent them. Both
% images and profiles are normalised to their own peak, which is what makes
% the shapes comparable.
for i = 1:4
    v = R{i}; m = max(abs(v(~isnan(v))));
    subplot(2,5,i+1);
    hI = imagesc(xax,xax,v,[-m m]); axis image off; colormap(gca,cm);
    set(hI,'AlphaData',~isnan(v)); set(gca,'Color',[0.42 0.47 0.52]);
    title(T{i},'FontSize',10);
    subplot(2,5,i+6);
    pr = v(31,:); pr = pr/max(abs(pr(~isnan(pr))));
    plot(xax,pr,'LineWidth',1.5); hold on; yline(0,'k:');
    axis square; grid on; ylim([-1.15 1.15]); yticks(-1:0.5:1);
    xlabel('x'); ylabel('normalised \Delta\sigma'); title('centre cut','FontSize',9);
end
exportgraphics(f,'fig_dagaf_separation.png','Resolution',150);
fprintf('wrote fig_dagaf_separation.png\n');
rmdir(tmpd,'s');

%% ---- figure 2: the DAGAF decomposition ------------------------------------
g = figure('Color','w','Position',[0 0 900 1000],'Visible','off');
K = size(di,1);
subplot(K+2,1,1); plot(X,'k','LineWidth',1); axis tight; grid on;
title('measured 208-voltage vector (16 injections \times 13 measurements)','FontSize',10);
set(gca,'XTick',0:13:208,'XTickLabel',[]);
for k = 1:K
    subplot(K+2,1,k+1); plot(di(k,:),'LineWidth',1); axis tight; grid on;
    ylabel(sprintf('IMF%d',k));
    if k==1, title('DAGAF components, nIMF=2 (window = 13 pts = one injection block)','FontSize',10); end
    set(gca,'XTick',0:13:208,'XTickLabel',[]);
end
subplot(K+2,1,K+2); plot(dr,'LineWidth',1); axis tight; grid on;
ylabel('res'); xlabel('measurement index'); set(gca,'XTick',0:13:208);
exportgraphics(g,'fig_dagaf_components.png','Resolution',150);
fprintf('wrote fig_dagaf_components.png\n');

function c = bwr(n)
h = floor(n/2); t = linspace(0,1,h)';
c = [[t t ones(h,1)]; [ones(h,1) flipud(t) flipud(t)]];
end
