%% EIT_FIG11_15_REPRO
% Reproduction of Figs. 11 and 15 of Cheng, Su & Ko, Curr Med Imaging 2022,
% with DAGAF added as a fourth panel.
%
%   Fig. 11  heart imaging: the lung-dominant component is REMOVED
%            (a) unprocessed, from zc_resp   -- it shows a ventilation cycle,
%                which zc_breathhold cannot contain
%            (b) PCA, PC3 removed | (c) EMD, IMF1 removed | (d) DAGAF, IMF1
%                removed -- all from zc_breathhold
%   Fig. 15  lung imaging: the lung-dominant component is BOOSTED, all panels
%            from zc_resp, matching the published figure where every panel
%            shows the same ventilation cycle.
%
% Settings follow PCA_montreal_data_1995_PCA_breathhld.m throughout: PCA on the
% 256-vector reshaped to 16x16, EMD/DAGAF on the 208 non-zero measurements with
% the zeros re-inserted, electrodes rotated [9:16,1:8], GREIT trained with
% lungs at 0.3, distr 3, Nsim 500, target_size 0.03, NF 0.5, meas_select on.

clc; close all;
PCREMOVE = 3;  PCNUM = 15;  BOOST = 1.01;   % 1% added, as in the 2022 paper

fmdl = mk_library_model('adult_male_16el_lungs');
fmdl.electrode = fmdl.electrode([9:16,1:8]);
img = mk_image(fmdl,1);
img.elem_data([fmdl.mat_idx{2}; fmdl.mat_idx{3}]) = 0.3;
[stim,msel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
img.fwd_model.stimulation = stim;
img.fwd_model = mdl_normalize(img.fwd_model,1);
imdl = mk_GREIT_model(img, 0.25, [], ...
   struct('imgsz',[32 32],'distr',3,'Nsim',500,'target_size',0.03, ...
          'target_offset',0,'noise_figure',0.5,'square_pixels',1));
imdl.fwd_model.meas_select = msel;

load montreal_data_1995
ZB = double(zc_breathhold);  ZR = double(zc_resp);
nzi = ZB(:,1)~=0;
gR = sum(ZR(nzi,:),1); [~,pkF] = max(abs(gR - gR(1)));
FR_R = max(2,pkF-8) : max(2,pkF-8)+15;          % 16 frames centred on the breath
FR_B = 2:17;
fprintf('zc_resp ventilation peak at frame %d -> frames %d-%d\n', pkF, FR_R(1), FR_R(end));

%% ---- process one recording with one weighting ---------------------------
function W = process(Z, nzi, method, kdrop, wval, PCREMOVE, PCNUM)
nF = size(Z,2);  W = zeros(256,nF);
for j = 1:nF
    vi = Z(:,j);
    switch method
      case 'PCA'
        X = zeros(16,16);
        for i=1:16, X(i,:) = vi((16*(i-1)+1):16*i).'; end
        Xm = mean(X,2);  X = X - Xm*ones(1,16);
        Y = X.'; Yc = Y - repmat(mean(Y,1),16,1);
        [~,~,Vv] = svd(Yc,'econ');  U = Vv;  PC = Yc*Vv;
        if wval == 0
            U(:,PCREMOVE) = [];  PC(:,PCREMOVE) = [];
            newX = U(:,1:PCNUM)*PC(:,1:PCNUM).';
        else
            PC(:,PCREMOVE) = PC(:,PCREMOVE)*wval;
            newX = U*PC.';
        end
        newX = newX + Xm*ones(1,16);
        out = vi;
        for i=1:16, out((16*(i-1)+1):16*i) = newX(i,:).'; end
      otherwise
        x = vi(nzi).';
        if strcmp(method,'EMD'), [imf,res] = emd_rcada(x);
        else, [imf,res] = dagaf2(x,2,2,'p','p',40,1e-8,3,6); end
        C = [imf;res];  w = ones(1,size(C,1));  w(kdrop) = wval;
        out = vi;  out(nzi) = (w*C).';
    end
    W(:,j) = out;
end
end

%% ---- render 4 panels of 16 frames ---------------------------------------
function render(imdl, DATA, REF, FRS, TTL, outfile)
tmp = tempname; mkdir(tmp); files = cell(1,4);
for k = 1:4
    im = inv_solve(imdl, REF{k}, DATA{k}(:,FRS{k}));
    im.calc_colours.ref_level = 0;  im.show_slices.img_cols = 4;
    h = figure('Color','w','Position',[0 0 540 560],'Visible','off');
    show_slices(im); axis off
    files{k} = fullfile(tmp,sprintf('p%d.png',k));
    exportgraphics(h, files{k},'Resolution',200,'ContentType','image'); close(h);
end
f = figure('Color','w','Position',[0 0 1720 520],'Visible','off');
for k = 1:4
    subplot(1,4,k); imshow(imread(files{k}));
    title(TTL{k},'FontSize',10,'FontWeight','bold','Interpreter','none');
    xlabel(sprintf('(%c)','a'+k-1),'FontSize',11);
end
exportgraphics(f, outfile,'Resolution',150); close(f); rmdir(tmp,'s');
fprintf('wrote %s\n', outfile);
end

%% ---- Fig. 11 : heart imaging -------------------------------------------
D11 = {ZR, process(ZB,nzi,'PCA',0,0,PCREMOVE,PCNUM), ...
           process(ZB,nzi,'EMD',1,0,PCREMOVE,PCNUM), ...
           process(ZB,nzi,'DAGAF',1,0,PCREMOVE,PCNUM)};
render(imdl, D11, {ZR(:,1),ZB(:,1),ZB(:,1),ZB(:,1)}, {FR_R,FR_B,FR_B,FR_B}, ...
   {'Original reconstruction (zc_resp)','PCA-based, PC3 removed', ...
    'EMD-based, IMF1 removed','DAGAF-based, IMF1 removed'}, 'fig11_repro_abcd.png');

%% ---- Fig. 15 : lung imaging --------------------------------------------
D15 = {ZR, process(ZR,nzi,'PCA',0,BOOST,PCREMOVE,PCNUM), ...
           process(ZR,nzi,'EMD',1,BOOST,PCREMOVE,PCNUM), ...
           process(ZR,nzi,'DAGAF',1,BOOST,PCREMOVE,PCNUM)};
render(imdl, D15, repmat({ZR(:,1)},1,4), repmat({FR_R},1,4), ...
   {'Original reconstruction', sprintf('PCA-based, %g%% PC3 added',(BOOST-1)*100), ...
    sprintf('EMD-based, %g%% IMF1 added',(BOOST-1)*100), ...
    sprintf('DAGAF-based, %g%% IMF1 added',(BOOST-1)*100)}, 'fig15_repro_abcd.png');
