function S = eit_human_separating_component()
%EIT_HUMAN_SEPARATING_COMPONENT  Same split as EIT_noise_geometry/src/
% p2_separating_component.m, applied to the human record of Section 3.7
% (montreal_data_1995, breath-hold, settings of Table 1 / eit_human_cmi_repro.m).
%
% GREIT one-step difference reconstruction with normalised data is affine in the
% measured frame for a fixed reference vref:  recon(v) = RM * ((v - vref)./vref).
% With removed = vi - contract_output,
%     recon(contract) = recon(vi) - RM*(removed./vref)
%                     =  [physiological term] + [artefact term]
% The artefact term is formed as  -inv_solve(imdl, vref, vref + removed)  and the
% closure of the split is checked against the directly reconstructed contract.
%
% Supplementary (as RESULTS_SEPARATING_COMPONENT.md section 5): the same contract
% applied to the reference frame alone, reconstructing what it removes there.
%
% No ground-truth heart ROI exists on human data. Statistics are reported over the
% whole image and over the "blob" region: pixels where |contract output| is at
% least half its maximum in that frame.
%
% Writes  human_separating_component.mat / .txt
%         fig_human_separating_component.png

clc; close all;
here = fileparts(mfilename('fullpath'));
outTxt = fullfile(here,'human_separating_component.txt');
if isfile(outTxt), delete(outTxt); end
diary(outTxt);
fprintf('eit_human_separating_component  %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));

PCREMOVE = 3; PCNUM = 15;      % as eit_human_cmi_repro.m (PC3 on human data)

%% --- model, exactly as eit_human_cmi_repro.m ------------------------------
fmdl = mk_library_model('adult_male_16el_lungs');
fmdl.electrode = fmdl.electrode([9:16,1:8]);
img = mk_image(fmdl, 1);
img.elem_data([fmdl.mat_idx{2}; fmdl.mat_idx{3}]) = 0.3;
[stim,msel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
img.fwd_model.stimulation = stim;
img.fwd_model = mdl_normalize(img.fwd_model, 1);
opt = struct('imgsz',[32 32],'distr',3,'Nsim',500,'target_size',0.03, ...
             'target_offset',0,'noise_figure',0.5,'square_pixels',1);
imdl = mk_GREIT_model(img, 0.25, [], opt);
imdl.fwd_model.meas_select = msel;

load montreal_data_1995
Zall = double(zc_breathhold);  nF = size(Zall,2);
nzi  = Zall(:,1) ~= 0;
vref = Zall(:,1);
fprintf('zc_breathhold: %d frames @ 7 fps, %d of 256 rows non-zero; reference = frame 1\n', nF, sum(nzi));

methods = {'PCA','EMD','DAGAF'};
frames  = 2:nF;
S = struct();
S.methods = methods; S.frames = frames;

%% --- contract output and removed part, per method and frame ---------------
for m = 1:3
    mth = methods{m};
    OUT = zeros(256,nF); REM = zeros(256,nF);
    for f = 1:nF
        vi = Zall(:,f);
        out = contract(vi, mth, nzi, PCREMOVE, PCNUM);
        OUT(:,f) = out; REM(:,f) = vi - out;
    end
    % what the same contract removes from the reference frame alone
    rem_ref = vref - contract(vref, mth, nzi, PCREMOVE, PCNUM);

    % measurement-domain sizes (non-zero rows only)
    rel_removed = sqrt(mean(sum(REM(nzi,frames).^2,1) ./ sum(Zall(nzi,frames).^2,1)));
    rel_phys    = sqrt(mean(sum((Zall(nzi,frames)-vref(nzi)).^2,1) ./ sum(Zall(nzi,frames).^2,1)));
    rel_remref  = norm(rem_ref(nzi))/norm(vref(nzi));
    % how much of removed(vi) is removed(vref): ||removed(vi)-removed(vref)||/||removed(vref)||
    rel_dyn     = sqrt(mean(sum((REM(nzi,frames)-rem_ref(nzi)).^2,1))) / norm(rem_ref(nzi));

    % reconstructions (frames x pixels)
    im_phys = inv_solve(imdl, vref, Zall(:,frames));           e_phys = im_phys.elem_data;   % unprocessed
    im_out  = inv_solve(imdl, vref, OUT(:,frames));            e_out  = im_out.elem_data;    % contract output
    im_rem  = inv_solve(imdl, vref, vref + REM(:,frames));     e_art  = -im_rem.elem_data;   % artefact term
    im_rref = inv_solve(imdl, vref, vref + rem_ref);           e_aref = -im_rref.elem_data;  % ref-frame part
    closure = norm(e_out - (e_phys + e_art), 'fro') / norm(e_out, 'fro');

    nfr = numel(frames);
    ratio_all = zeros(nfr,1); ratio_blob = zeros(nfr,1); r_art = zeros(nfr,1); r_phys = zeros(nfr,1);
    ref_share = zeros(nfr,1); r_aref = zeros(nfr,1); fom_art_blob = zeros(nfr,1); nblob = zeros(nfr,1);
    for k = 1:nfr
        o = e_out(:,k); p = e_phys(:,k); a = e_art(:,k);
        blob = abs(o) >= 0.5*max(abs(o));  nblob(k) = sum(blob);
        ratio_all(k)  = rms(a)/rms(p);
        ratio_blob(k) = rms(a(blob))/rms(p(blob));
        r_art(k)  = pcorr(o, a);  r_phys(k) = pcorr(o, p);
        ref_share(k) = rms(e_aref(blob))/rms(a(blob));
        r_aref(k) = pcorr(a, e_aref);
    end
    S.(mth) = struct('rel_removed',rel_removed,'rel_phys',rel_phys,'rel_remref',rel_remref, ...
        'rel_dyn',rel_dyn,'closure',closure,'ratio_all',ratio_all,'ratio_blob',ratio_blob, ...
        'r_art',r_art,'r_phys',r_phys,'ref_share',ref_share,'r_aref',r_aref,'nblob',nblob, ...
        'e_phys',e_phys,'e_out',e_out,'e_art',e_art,'e_aref',e_aref);

    fprintf('\n=== %s (removes %s) ===\n', mth, idxname(mth,PCREMOVE));
    fprintf('measurement domain: ||removed||/||vi||  %.3f (rms over frames)   ||vi-vref||/||vi||  %.4f   ratio %.1f\n', ...
            rel_removed, rel_phys, rel_removed/rel_phys);
    fprintf('                    ||removed(vref)||/||vref|| %.3f ;  ||removed(vi)-removed(vref)||/||removed(vref)|| %.3f\n', rel_remref, rel_dyn);
    fprintf('split closure ||out-(phys+art)||/||out|| = %.2e\n', closure);
    fprintf('image, %d frames:  ratio_all  rms(art)/rms(phys)   median %.1f  [min %.1f, max %.1f]\n', nfr, median(ratio_all), min(ratio_all), max(ratio_all));
    fprintf('                   ratio_blob (|out|>=0.5 max)     median %.1f  [min %.1f, max %.1f]   blob size median %d px\n', median(ratio_blob), min(ratio_blob), max(ratio_blob), median(nblob));
    fprintf('                   r(out, art)  median %.3f [%.3f, %.3f];  r(out, phys) median %.3f [%.3f, %.3f]\n', ...
            median(r_art), min(r_art), max(r_art), median(r_phys), min(r_phys), max(r_phys));
    fprintf('reference-frame part: rms(aref)/rms(art) in blob  median %.2f [%.2f, %.2f];  r(art, aref) median %.3f [%.3f, %.3f]\n', ...
            median(ref_share), min(ref_share), max(ref_share), median(r_aref), min(r_aref), max(r_aref));
end

%% --- figure: frame 10 (one of the published frames), one scale per row ----
fshow = 10; k = find(frames==fshow);
fig = figure('Color','w','Position',[0 0 1400 900],'Visible','off');
cols = {'unprocessed (= physiological term)','contract output','-recon(removed): artefact term','-recon(removed from vref alone)'};
for m = 1:3
    R = S.(methods{m});
    E = {R.e_phys(:,k), R.e_out(:,k), R.e_art(:,k), R.e_aref};
    cl = max(abs(R.e_out(:,k)));
    for c = 1:4
        subplot(3,4,(m-1)*4+c);
        tmp = im_out; tmp.elem_data = E{c};
        sl = calc_slices(tmp); 
        imagesc(sl, [-cl cl]); axis image off; colormap(bwr(256));
        if m==1, title(cols{c},'FontSize',9); end
        if c==1, text(-0.15,0.5,methods{m},'Units','normalized','Rotation',90,'HorizontalAlignment','center','FontWeight','bold'); end
        text(0.02,0.04,sprintf('rms %.3g',rms(E{c})),'Units','normalized','FontSize',7);
    end
end
sgtitle(sprintf('montreal_data_1995 breath-hold, frame %d; each row on the scale of its own contract output', fshow),'Interpreter','none');
exportgraphics(fig, fullfile(here,'fig_human_separating_component.png'), 'Resolution',150);
fprintf('\nwrote fig_human_separating_component.png (frame %d)\n', fshow);
save(fullfile(here,'human_separating_component.mat'), '-struct', 'S');
diary off
end

%% ------------------------------------------------------------------------
function out = contract(vi, mth, nzi, PCREMOVE, PCNUM)
switch mth
    case 'PCA'      % 16x16 on the full 256-vector, as the published human pipeline
        X = zeros(16,16);
        for i = 1:16, X(i,:) = vi((16*(i-1)+1):16*i).'; end
        Xo_mean = mean(X,2);
        X = X - Xo_mean*ones(1,16);
        [U,PC] = princomp_(X.');
        U(:,PCREMOVE) = [];  PC(:,PCREMOVE) = [];
        newX = U(:,1:PCNUM)*PC(:,1:PCNUM).' + Xo_mean*ones(1,16);
        out = zeros(256,1);
        for i = 1:16, out((16*(i-1)+1):16*i) = newX(i,:).'; end
    otherwise       % EMD / DAGAF on the 208 non-zero rows, drop mode 1
        x = vi(nzi).';
        if strcmp(mth,'EMD'), [imf,res] = emd_rcada(x);
        else,                 [imf,res] = dagaf2(x,2,2,'p','p',40,1e-8,3,6); end
        C = [imf; res];  w = ones(1,size(C,1)); w(1) = 0;
        out = vi;  out(nzi) = (w*C).';
end
end

function s = idxname(mth, PCREMOVE)
if strcmp(mth,'PCA'), s = sprintf('PC%d of 16x16 reshape', PCREMOVE); else, s = 'mode 1'; end
end

function [COEFF,SCORE,latent] = princomp_(Y)
mu = mean(Y,1);  Yc = Y - repmat(mu,size(Y,1),1);
[~,Sg,V] = svd(Yc,'econ');
COEFF = V;  SCORE = Yc*V;  latent = diag(Sg).^2/max(size(Y,1)-1,1);
end

function r = pcorr(a, b)
R = corrcoef(a(:), b(:)); r = R(1,2);
end

function c = bwr(n)
h=floor(n/2); q=linspace(0,1,h)';
c=[[q q ones(h,1)];[ones(h,1) flipud(q) flipud(q)]];
end
