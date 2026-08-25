%% EIT_GRAPHICAL_ABSTRACT
% One image carrying the paper's argument: a single frame of 208 voltages is
% decomposed under a contract that every method must satisfy, one row is
% removed, and the same reconstruction follows. The three methods then differ
% only in the decomposition, and ringing is what separates them.
%
% Writes graphical_abstract.png.

clc; clear; close all;

P.r_lung = 0.6;  P.r_heart = 0.3;  P.sig_lung = 0.5;  P.sig_heart = 2.0;
P.pca_keep = 1:12;  P.pca_key = 4;  P.emd_key = 1;
P.dagaf = struct('nIMF',6,'chi',2,'alpha',3,'eType','p','maskCap',6,'TH1',40,'TH2',1e-8);

here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(here,'rcada_emd')); addpath(fullfile(here,'github_release','src'));
stim = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);

fprintf('phantom and GREIT model ...\n');
lungs = sprintf('solid obj1 = sphere(1,0.1,1;%.2f) or sphere(-1,0.1,1;%.2f);', P.r_lung, P.r_lung);
heart = sprintf('solid obj2 = sphere(0,0.1,1;%.2f);', P.r_heart);
[fmdl, midx] = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1],{'obj1','obj2',[lungs heart]});
fmdl.stimulation = stim;
img = mk_image(fmdl,1);  vh = fwd_solve(img);
img.elem_data(midx{2}) = P.sig_lung;  img.elem_data(midx{3}) = P.sig_heart;
vi = fwd_solve(img);  X = double(vi.meas(:)).';

opt = struct('imgsz',[64 64],'distr',0,'Nsim',1000,'target_size',0.02,'noise_figure',0.5);
fg = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1]);  fg.stimulation = stim;  fg = mdl_normalize(fg,1);
imdl = mk_GREIT_model(mk_image(fg,1),0.2,[],opt);

% --- decompositions ------------------------------------------------------
A = reshape(X,13,16); mu = mean(A,1); Ac = A - repmat(mu,13,1);
[V,D] = eig(cov(Ac)); [~,o] = sort(diag(D),'descend'); V = V(:,o);
PC = Ac*V; K = size(V,2);
Cp = zeros(K+1,208);
for k = 1:K, Cp(k,:) = reshape(PC(:,k)*V(:,k)',1,[]); end
Cp(K+1,:) = reshape(repmat(mu,13,1),1,[]);
[imf,res] = emd_rcada(X);  Ce = [imf; res];
[i2,r2]  = dagaf2(X,P.dagaf.nIMF,P.dagaf.chi,P.dagaf.eType,P.dagaf.eType, ...
                  P.dagaf.TH1,P.dagaf.TH2,P.dagaf.alpha,P.dagaf.maskCap);
Cd = [i2; r2];

keep = false(1,size(Cp,1)); keep([P.pca_keep size(Cp,1)]) = true;
wp = double(keep);            wp(P.pca_key) = 0;
we = ones(1,size(Ce,1));      we(P.emd_key) = 0;
wd = ones(1,size(Cd,1));      wd(1)         = 0;

R = @(C,w) recon(imdl,vh,vi,C,w);
I = {inv_solve(imdl,vh,vi), R(Cp,wp), R(Ce,we), R(Cd,wd)};
for j = 1:4, I{j}.calc_colours.ref_level = 0; end
tgt = [0;0.1;1;P.r_heart];
rng_ = zeros(1,4);
for j = 1:4, p = eval_GREIT_fig_merit(I{j},tgt); rng_(j) = p(5); end
fprintf('  RNG: %.2f  %.3f  %.3f  %.3f\n', rng_);

%% ---- layout: 16:9, following the journal's own graphical-abstract style --
f = figure('Color','w','Position',[0 0 1333 750],'Visible','off');
box_ = @(x,y,w,h,str,fs) annotation('textbox',[x y w h],'String',str, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',fs,'LineWidth',1,'EdgeColor',[.2 .2 .2],'BackgroundColor','w','Margin',3);
arr = @(x1,y1,x2,y2) annotation('arrow',[x1 x2],[y1 y2],'LineWidth',1.4,'Color',[.25 .40 .70]);

% ---- row 1: phantom -> forward problem -> measured frame ----------------
if exist('fig1a_phantom.png','file')
    axes('Position',[0.025 0.585 0.235 0.355]);
    imshow(imread('fig1a_phantom.png')); axis off
end
text(0.5,1.06,'3-D FEM phantom','Units','normalized','HorizontalAlignment','center', ...
     'FontSize',11,'FontWeight','bold');

box_(0.272,0.715,0.105,0.085,'EIT forward problem',9);
arr(0.265,0.757,0.270,0.757);  arr(0.379,0.757,0.398,0.757);

ax = axes('Position',[0.405 0.615 0.265 0.295]);
plot(1:208, X, 'k-','LineWidth',0.8); hold on
for b = 13:13:195, xline(b,':','Color',[.78 .78 .78]); end
xlim([1 208]); box on; set(ax,'FontSize',8,'XTick',[1 104 208]);
xlabel('measurement index','FontSize',8);
title('One frame: 208 measured voltages','FontSize',11,'FontWeight','bold');

% ---- the contract, fed to three decompositions --------------------------
box_(0.700,0.630,0.275,0.265, ...
    {'\bfThe common contract\rm','', ...
     'each decomposition returns rows', 'summing exactly to the frame', ...
     '\Sigma_k C_k = x', '', 'separation removes one row:  w_j = 0'}, 10);
arr(0.678,0.762,0.698,0.762);

names3 = {'PCA','EMD','DAGAF'};
for j = 1:3
    box_(0.715+(j-1)*0.090, 0.500, 0.078, 0.062, ['\bf' names3{j}], 10);
    arr(0.754+(j-1)*0.090, 0.625, 0.754+(j-1)*0.090, 0.566);
end

% ---- row 2: the four reconstructions ------------------------------------
annotation('textbox',[0.03 0.500 0.62 0.055],'String', ...
   '\bfHeart imaging from a single frame\rm  (the heart is invisible before processing)', ...
   'EdgeColor','none','HorizontalAlignment','left','FontSize',11);
name = {'Unprocessed','PCA-based','EMD-based','DAGAF-based'};
x0 = 0.045; w = 0.135; gap = 0.020;
for j = 1:4
    axes('Position',[x0+(j-1)*(w+gap) 0.215 w 0.275]);
    show_slices(I{j}); axis image off
    title(name{j},'FontSize',10.5,'FontWeight','bold');
    c = [.25 .25 .25]; fw = 'normal';
    if j == 4, c = [.75 .10 .10]; fw = 'bold'; end
    text(0.5,-0.11, sprintf('ringing  %.2f', rng_(j)),'Units','normalized', ...
        'HorizontalAlignment','center','FontSize',10,'Color',c,'FontWeight',fw);
end

% ---- findings -----------------------------------------------------------
annotation('textbox',[0.685 0.205 0.290 0.290],'String', ...
  {'\bfFindings:\rm', ...
   '\bullet  A common contract makes measurement-domain', ...
   '    decompositions comparable on equal terms.', ...
   '\bullet  Ringing separates the three methods; shape', ...
   '    deformation cannot, being the selection criterion.', ...
   '\bullet  A fixed convolution with no adaptive element', ...
   '    reproduces DAGAF to within 1.2%.', ...
   '\bullet  The removed component is not the lung component.'}, ...
  'EdgeColor',[.75 .75 .75],'BackgroundColor',[.97 .97 .98], ...
  'VerticalAlignment','top','FontSize',9.5,'Margin',6);

% ---- takeaway strip -----------------------------------------------------
annotation('textbox',[0.025 0.038 0.950 0.115],'String', ...
  ['\bfWhat distinguishes the three methods is whether the removed component is ' ...
   're-derived from the frame\rm - not the sophistication of the decomposition.'], ...
  'EdgeColor',[.25 .40 .70],'BackgroundColor',[.93 .95 .99],'LineWidth',1.2, ...
  'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);

exportgraphics(f,'graphical_abstract.png','Resolution',220);
fprintf('wrote graphical_abstract.png\n');

%% ========================================================================
function img = recon(imdl, vh, vi, comp, w)
v = vi; v.meas = (w(:).'*comp).';
img = inv_solve(imdl, vh, v); img.calc_colours.ref_level = 0;
end
function o = ternary(c,a,b), if c, o=a; else, o=b; end, end
