%% EIT_PARTIAL_WEIGHT
% The CMI 2022 method removes the lung-dominant component outright (w = 0)
% for heart imaging. That is a binary choice; nothing requires it. Here the
% weight on that component is swept continuously from 1 (no processing) to 0
% (full removal), for all three decompositions.
%
% Motivation: DAGAF matches PCA on shape (SD 0.030 vs 0.025) but is far worse
% on ringing (RNG 7.77 vs 1.46). If the ringing comes from removing too much
% of IMF1 -- which carries the sharp peaks at the injection-block boundaries,
% and hence the heart's high spatial frequencies -- then a partial weight
% should trade SD against RNG and may find a better operating point.
%
% Requires startup.m (Netgen + EIDORS + rcada_emd on the path).

clc; close all;
P.r_lung = 0.6; P.r_heart = 0.3; P.sig_lung = 0.5; P.sig_heart = 2.0;
W = 0:0.05:1;

stim = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
lungs = sprintf('solid obj1 = sphere(1,0.1,1;%.2f) or sphere(-1,0.1,1;%.2f);',P.r_lung,P.r_lung);
heart = sprintf('solid obj2 = sphere(0,0.1,1;%.2f);',P.r_heart);
[fmdl,midx] = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1],{'obj1','obj2',[lungs heart]});
fmdl.stimulation = stim;
im = mk_image(fmdl,1); vh = fwd_solve(im);
im.elem_data(midx{2}) = P.sig_lung; im.elem_data(midx{3}) = P.sig_heart;
vi = fwd_solve(im);
X = double(vi.meas(:)).';

opt = struct('imgsz',[64 64],'distr',0,'Nsim',1000,'target_size',0.02,'noise_figure',0.5);
fg = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1]); fg.stimulation = stim; fg = mdl_normalize(fg,1);
imdl = mk_GREIT_model(mk_image(fg,1),0.2,[],opt);
tgt = [0;0.1;1;P.r_heart];

%% --- the four configurations -------------------------------------------
A = reshape(X(:),13,16); mu = mean(A,1); Ac = A-repmat(mu,13,1);
[V,L] = eig(cov(Ac)); [~,o] = sort(diag(L),'descend'); V = V(:,o); PCs = Ac*V;
Cp = zeros(17,208);
for k=1:16, Cp(k,:) = reshape(PCs(:,k)*V(:,k)',1,[]); end
Cp(17,:) = reshape(repmat(mu,13,1),1,[]);

[ei,er] = emd_rcada(X);
[d6i,d6r] = dagaf2(X,2,2,'p','p',40,1e-8,3,6);
[d4i,d4r] = dagaf2(X,2,2,'p','p',40,1e-8,3,4);

CFG = { struct('name','PCA (PC4)',        'C',Cp,        'k',4)
        struct('name','EMD (IMF1)',       'C',[ei;er],   'k',1)
        struct('name','DAGAF cap6 (IMF1)','C',[d6i;d6r], 'k',1)
        struct('name','DAGAF cap4 (IMF1)','C',[d4i;d4r], 'k',1) };

%% --- sweep ---------------------------------------------------------------
R = nan(numel(CFG), numel(W), 3);
for c = 1:numel(CFG)
    C = CFG{c}.C;  K = size(C,1);  kk = CFG{c}.k;
    fprintf('\n%s\n%6s %8s %8s %9s\n', CFG{c}.name, 'w', 'PE','SD','RNG');
    for j = 1:numel(W)
        w = ones(1,K); w(kk) = W(j);
        v = vi; v.meas = (w*C).';
        p = eval_GREIT_fig_merit(inv_solve(imdl,vh,v), tgt);
        R(c,j,:) = p([2 4 5]);
        fprintf('%6.2f %8.3f %8.3f %9.3f\n', W(j), p(2), p(4), p(5));
    end
end

%% --- best operating point per method ------------------------------------
fprintf('\n================ best operating point (SD < 0.10) ================\n');
fprintf('%-18s %6s %8s %8s %9s\n','method','w','PE','SD','RNG');
for c = 1:numel(CFG)
    sd = squeeze(R(c,:,2)); rng_ = squeeze(R(c,:,3));
    ok = find(sd < 0.10);
    if isempty(ok), fprintf('%-18s   (no w with SD < 0.10)\n', CFG{c}.name); continue; end
    [~,ii] = min(rng_(ok)); j = ok(ii);
    fprintf('%-18s %6.2f %8.3f %8.3f %9.3f\n', CFG{c}.name, W(j), R(c,j,1), R(c,j,2), R(c,j,3));
    fprintf('%-18s %6.2f %8.3f %8.3f %9.3f   <- binary removal (w=0)\n','   (vs)',0,R(c,1,1),R(c,1,2),R(c,1,3));
end

%% --- figure --------------------------------------------------------------
% SD spans 0.03 to 0.95 across the sweep, and the values that matter -- the
% four methods at w = 0 -- sit between 0.030 and 0.059. On a linear axis they
% collapse into the bottom few per cent, so SD gets the same log treatment as
% RNG, in both the sweep and the trade-off panel.
f = figure('Color','w','Position',[0 0 1320 470],'Visible','off');
co = lines(4);
subplot(1,3,1); hold on
for c=1:4
    plot(W,squeeze(R(c,:,2)),'LineWidth',1.8,'Color',co(c,:));
    plot(0,R(c,1,2),'o','MarkerSize',6,'MarkerFaceColor',co(c,:),'MarkerEdgeColor','k');
end
xlabel('weight on the lung-dominant component'); ylabel('SD'); grid on
set(gca,'YScale','log'); ylim([0.01 1.2]); axis square
title({'shape deformation (log scale)','markers = binary removal'});

subplot(1,3,2); hold on
for c=1:4, plot(W,squeeze(R(c,:,3)),'LineWidth',1.8,'Color',co(c,:)); end
xlabel('weight'); ylabel('RNG'); grid on; set(gca,'YScale','log')
title('ringing (log scale)'); axis square
legend(cellfun(@(s)s.name,CFG,'uni',0),'Location','northwest'); legend boxoff

subplot(1,3,3); hold on
for c=1:4
    plot(squeeze(R(c,:,2)),squeeze(R(c,:,3)),'-','LineWidth',1.5,'Color',co(c,:));
    plot(R(c,1,2),R(c,1,3),'p','MarkerSize',15,'MarkerFaceColor',co(c,:), ...
         'MarkerEdgeColor','k');
    text(R(c,1,2)*1.10, R(c,1,3), sprintf('%.3f',R(c,1,2)), 'FontSize',8.5);
end
xlabel('SD'); ylabel('RNG'); grid on
set(gca,'XScale','log','YScale','log'); xlim([0.02 1.2]); ylim([0.8 60]);
title({'trade-off curve','star = binary removal, w = 0'}); axis square
exportgraphics(f,'fig_partial_weight.png','Resolution',150);
fprintf('\nwrote fig_partial_weight.png\n');
save('partial_weight_sweep.mat','R','W','CFG');
