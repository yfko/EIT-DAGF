%% EIT_LUNG_TABLE
% Lung-imaging GREIT figures of merit, reproduced against the manuscript
% "Error quantification_GREIT Test Parameters" and then recomputed with a
% mask that actually isolates one lung.
%
% Three evaluation variants per method:
%   (A) manuscript, no negation   : elem_data(1614:end)=0
%   (B) manuscript, with negation : elem_data(1614:end)=0, then negate
%   (C) corrected                 : negate, then keep x > 0 (one lung)
%
% (B) is what the manuscript settled on. It fixes calc_ringing's img<0
% assumption but keeps a mask that retains the BOTTOM HALF in y -- both
% lungs' lower halves plus the heart's -- rather than one lung.

clc; close all;
P.r_heart = 0.3; P.sig_lung = 0.5; P.sig_heart = 2.0;
P.w_boost = 1.001;            % as in the manuscript's phantom script
RL = [0.3 0.6];

% manuscript values, for side-by-side printing: {r_lung}{variant}(metric, method)
MS.r03.B = [0.75 0.75 0.76; 0.96 0.96 0.93; 0.45 0.44 0.38];
MS.r03.A = [0.69 0.37 0.55; 0.12 0.15 0.16; 3.91 3.37 4.40];
MS.r06.B = [0.74 0.74 0.74; 0.72 0.71 0.69; 0.06 0.07 0.07];
MS.r06.A = [0.04 0.06 0.04; 0.14 0.13 0.17; 29.49 23.75 27.03];

stim = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
opt  = struct('imgsz',[64 64],'distr',0,'Nsim',1000,'target_size',0.02,'noise_figure',0.5);
fg = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1]); fg.stimulation = stim; fg = mdl_normalize(fg,1);
imdl = mk_GREIT_model(mk_image(fg,1),0.2,[],opt);
ct  = interp_mesh(imdl.rec_model);
pix = (ct(1:2:end,:) + ct(2:2:end,:))/2;

for ri = 1:numel(RL)
    rl = RL(ri);
    lungs = sprintf('solid obj1 = sphere(1,0.1,1;%.2f) or sphere(-1,0.1,1;%.2f);',rl,rl);
    heart = sprintf('solid obj2 = sphere(0,0.1,1;%.2f);',P.r_heart);
    [fm,midx] = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1],{'obj1','obj2',[lungs heart]});
    fm.stimulation = stim;
    im = mk_image(fm,1); vh = fwd_solve(im);
    im.elem_data(midx{2}) = P.sig_lung; im.elem_data(midx{3}) = P.sig_heart;
    vi = fwd_solve(im);
    X = double(vi.meas(:)).';
    tgt = [1; 0.1; 1; rl];

    % --- decompositions, lung-boosted recomposition ----------------------
    A = reshape(X(:),13,16); mu = mean(A,1); Ac = A - repmat(mu,13,1);
    [V,L] = eig(cov(Ac)); [~,o] = sort(diag(L),'descend'); V = V(:,o); PCs = Ac*V;
    Cp = zeros(17,208);
    for k=1:16, Cp(k,:) = reshape(PCs(:,k)*V(:,k)',1,[]); end
    Cp(17,:) = reshape(repmat(mu,13,1),1,[]);
    [ei,er] = emd_rcada(X);
    [di,dr] = dagaf2(X,2,2,'p','p',40,1e-8,3,6);

    IMGS = {inv_solve(imdl,vh,vi), ...
            boost(imdl,vh,vi,Cp,      4, P.w_boost), ...
            boost(imdl,vh,vi,[ei;er], 1, P.w_boost), ...
            boost(imdl,vh,vi,[di;dr], 1, P.w_boost)};
    NAMES = {'Original','PCA-based','EMD-based','DAGAF-based'};

    fprintf('\n\n#########  LUNG IMAGING,  r_lung = %.1f  #########\n', rl);
    ms = MS.(sprintf('r%02d',round(rl*10)));
    variants = {'A: manuscript mask, NO negation', ...
                'B: manuscript mask + negation  ', ...
                'C: corrected (negate, keep x>0)'};
    for v = 1:3
        fprintf('\n--- %s ---\n', variants{v});
        fprintf('%-6s %11s %11s %11s %13s\n','', NAMES{1},NAMES{2},NAMES{3},NAMES{4});
        rows = {'PE','SD','RNG'};
        Vals = nan(3,4);
        for m = 1:4
            p = fom(IMGS{m}, v, pix, tgt);
            Vals(:,m) = p([2 4 5]);
        end
        for i = 1:3
            fprintf('%-6s', rows{i});
            fprintf('%11.3f%11.3f%11.3f%13.3f\n', Vals(i,1),Vals(i,2),Vals(i,3),Vals(i,4));
            if v <= 2
                f = 'A'; if v==2, f = 'B'; end
                fprintf('%-6s%11.2f%11.2f%11.2f%13s   <- manuscript\n', ...
                        '', ms.(f)(i,1), ms.(f)(i,2), ms.(f)(i,3), '-');
            end
        end
    end
end

%% ------------------------------------------------------------------------
function img = boost(imdl, vh, vi, C, k, w)
ww = ones(1,size(C,1)); ww(k) = w;
v = vi; v.meas = (ww*C).';
img = inv_solve(imdl, vh, v); img.calc_colours.ref_level = 0;
end
function p = fom(img, variant, pix, tgt)
a = img;
switch variant
    case 1, a.elem_data(3228/2:end) = 0;
    case 2, a.elem_data(3228/2:end) = 0; a.elem_data = -a.elem_data;
    case 3, a.elem_data = -a.elem_data; a.elem_data(pix(:,1) <= 0) = 0;
end
p = eval_GREIT_fig_merit(a, tgt);
end
