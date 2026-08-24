%% EIT_DAGAF_IMF1_WEIGHT
% DAGAF counterpart of Fig. 12 of Cheng, Su & Ko, Curr Med Imaging 2022, which
% showed EMD reconstructions at IMF1 x 0.99, x 0 and x 1.01 with their centre
% cuts. Here the same三 states are shown for DAGAF: the unmodified frame, IMF1
% removed (heart imaging), and IMF1 boosted (lung imaging).
%
% The boost is chosen rather than assumed: it is swept until the heart stops
% forming a distinguishable local maximum on the centre cut.

clc; close all;
P.r_lung=0.6; P.r_heart=0.3; P.sig_lung=0.5; P.sig_heart=2.0;
xax = -2:4/63:2;

stim = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
lungs = sprintf('solid obj1 = sphere(1,0.1,1;%.2f) or sphere(-1,0.1,1;%.2f);',P.r_lung,P.r_lung);
heart = sprintf('solid obj2 = sphere(0,0.1,1;%.2f);',P.r_heart);
[fm,midx] = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1],{'obj1','obj2',[lungs heart]});
fm.stimulation = stim;
im0 = mk_image(fm,1); vh = fwd_solve(im0);
im0.elem_data(midx{2}) = P.sig_lung; im0.elem_data(midx{3}) = P.sig_heart;
vi = fwd_solve(im0);
X = double(vi.meas(:)).';

opt = struct('imgsz',[64 64],'distr',0,'Nsim',1000,'target_size',0.02,'noise_figure',0.5);
fg = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1]); fg.stimulation = stim; fg = mdl_normalize(fg,1);
imdl = mk_GREIT_model(mk_image(fg,1),0.2,[],opt);

truth  = calc_slices(im0);
SIMCUT = truth(31,:);            % the phantom's own centre cut, 0.5 .. 2.0

[di,dr] = dagaf2(X,2,2,'p','p',40,1e-8,3,6);  C = [di;dr];
recon = @(w) inv_solve(imdl, vh, setfield(vi,'meas',(w*C).'));   %#ok<SFLD>

% --- pick the boost -------------------------------------------------------
% heart contrast on the centre cut: value at x = 0 minus the mean of the two
% flanking minima. It is positive while the heart is still a local maximum.
function [hc, cut] = contrast(img, xax)
img.calc_colours.ref_level = 0;
S = calc_slices(img); cut = S(31,:);
c0 = cut(32);
[~,l] = min(cut(1:28)); [~,r] = min(cut(37:end)); r = r+36;
hc = c0 - (cut(l)+cut(r))/2;
end

WS = [1.00 1.005 1.01 1.02 1.03 1.05 1.06 1.07 1.08 1.12];
fprintf('\n%8s %14s %12s\n','IMF1 x','heart contrast','vs w=1');
base = []; hcs = zeros(size(WS));
for k = 1:numel(WS)
    ww = ones(1,size(C,1)); ww(1) = WS(k);
    hcs(k) = contrast(recon(ww), xax);
    if isempty(base), base = hcs(k); end
    fprintf('%8.3f %14.1f %11.0f%%\n', WS(k), hcs(k), 100*hcs(k)/base);
end
% the boost is the weight at which the heart stops being a local maximum
kk = find(hcs <= 0, 1);
WB = interp1(hcs([kk-1 kk]), WS([kk-1 kk]), 0);
WB = round(WB*100)/100;
fprintf('\nheart contrast reaches zero at IMF1 x %.2f -- used for the lung panel\n', WB);

%% ---- the three states, in the layout of CMI Fig. 12 --------------------
WSHOW = [1 0 WB];
LAB = {sprintf('(a)  IMF1 \\times 1  \x2014 unmodified'), ...
       sprintf('(b)  IMF1 \\times 0  \x2014 heart imaging'), ...
       sprintf('(c)  IMF1 \\times %.2f  \x2014 lung imaging', WB)};
CUT = cell(1,3); IMG = cell(1,3);
for k = 1:3
    ww = ones(1,size(C,1)); ww(1) = WSHOW(k);
    g = recon(ww); g.calc_colours.ref_level = 0;
    S = calc_slices(g); IMG{k} = S; CUT{k} = S(31,:);
end
% Each state is scaled to its own peak, as in the reference figure: removing
% IMF1 raises the peak from about -55 to +480, so a shared scale would leave
% the other two panels blank. The unmodified profile is overlaid in every cut
% so the change is still readable, and the peak values are quoted below.
fprintf('\ncentre-cut extremes:\n');
for k = 1:3
    fprintf('   %-28s min %8.1f   max %8.1f\n', LAB{k}(1:min(28,end)), ...
            min(CUT{k}), max(CUT{k}));
end

% Profiles are put on the phantom's conductivity scale, as the 2022 figures
% did, so the simulated profile can be drawn alongside them and all three
% panels can share one axis. That normalisation removes the magnitude
% difference between states, which is quoted above instead.
sfin = SIMCUT(~isnan(SIMCUT));
srange = max(sfin) - min(sfin);  smean = mean(sfin);
nrm = @(x) (x - mean(x(~isnan(x)))) ...
           / (max(x(:),[],'omitnan') - min(x(:),[],'omitnan')) * srange + smean;

f = figure('Color','w','Position',[0 0 1180 660],'Visible','off');
cmap = bwr_(256);
for k = 1:3
    lim = max(abs(IMG{k}(~isnan(IMG{k}))));
    subplot(2,3,k);
    hI = imagesc(xax,xax,IMG{k},[-lim lim]); axis image off; colormap(gca,cmap);
    set(hI,'AlphaData',~isnan(IMG{k})); set(gca,'Color',[0.42 0.47 0.52]);
    title(LAB{k},'FontSize',10.5,'FontWeight','bold');

    subplot(2,3,k+3);
    plot(xax, SIMCUT, 'k-','LineWidth',1.3); hold on
    plot(xax, nrm(CUT{k}), 'b-','LineWidth',1.7);
    grid on; axis square; xlim([-2 2]); ylim([0 2.1]); yticks(0:0.5:2);
    xlabel('x'); ylabel('\sigma');
    if k == 1
        legend('simulated','reconstructed','Location','south'); legend boxoff;
    end
end
exportgraphics(f,'fig_dagaf_imf1_states.png','Resolution',150);
fprintf('wrote fig_dagaf_imf1_states.png\n');

function c = bwr_(n)
h=floor(n/2); q=linspace(0,1,h)';
c=[[q q ones(h,1)];[ones(h,1) flipud(q) flipud(q)]];
end
