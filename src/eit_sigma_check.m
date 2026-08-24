%% EIT_SIGMA_CHECK
% Which heart conductivity produced the published Table 1?
%   "Error quantification_GREIT Test Parameters.m" line 52 uses 1.5
%   ("暫定為2改為1.5使用"), while thorax_GREIT.m line 17 and the CMI 2022
%   text both say 2. SD and RNG reproduce well either way; PE is the
%   sensitive one, so it should discriminate.
%
% Published Table 1 (heart imaging, r_lung = 0.6):
%   PE  0.04 / 0.04 / 0.04     SD 0.98 / 0.03 / 0.02     RNG 38.77 / 1.45 / 0.77
%   (Original / PCA-based / EMD-based)

clc; close all;
P.r_lung = 0.6; P.r_heart = 0.3; P.sig_lung = 0.5;
SIG_H = [1.5 2.0];
PUB = [0.04 0.04 0.04; 0.98 0.03 0.02; 38.77 1.45 0.77];

stim = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
opt  = struct('imgsz',[64 64],'distr',0,'Nsim',1000,'target_size',0.02,'noise_figure',0.5);
fg = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1]); fg.stimulation = stim; fg = mdl_normalize(fg,1);
imdl = mk_GREIT_model(mk_image(fg,1),0.2,[],opt);
tgt = [0; 0.1; 1; P.r_heart];

lungs = sprintf('solid obj1 = sphere(1,0.1,1;%.2f) or sphere(-1,0.1,1;%.2f);',P.r_lung,P.r_lung);
heart = sprintf('solid obj2 = sphere(0,0.1,1;%.2f);',P.r_heart);
[fm,midx] = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1],{'obj1','obj2',[lungs heart]});
fm.stimulation = stim;

for si = 1:numel(SIG_H)
    sh = SIG_H(si);
    im = mk_image(fm,1); vh = fwd_solve(im);
    im.elem_data(midx{2}) = P.sig_lung; im.elem_data(midx{3}) = sh;
    vi = fwd_solve(im);
    X = double(vi.meas(:)).';

    A = reshape(X(:),13,16); mu = mean(A,1); Ac = A - repmat(mu,13,1);
    [V,L] = eig(cov(Ac)); [~,o] = sort(diag(L),'descend'); V = V(:,o); PCs = Ac*V;
    Cp = zeros(17,208);
    for k=1:16, Cp(k,:) = reshape(PCs(:,k)*V(:,k)',1,[]); end
    Cp(17,:) = reshape(repmat(mu,13,1),1,[]);
    [ei,er] = emd_rcada(X);
    [di,dr] = dagaf2(X,2,2,'p','p',40,1e-8,3,6);

    IM = {inv_solve(imdl,vh,vi), drop(imdl,vh,vi,Cp,4), ...
          drop(imdl,vh,vi,[ei;er],1), drop(imdl,vh,vi,[di;dr],1)};

    V4 = nan(3,4);
    for m = 1:4
        p = eval_GREIT_fig_merit(IM{m}, tgt);
        V4(:,m) = p([2 4 5]);
    end

    fprintf('\n=========  sigma_heart = %.1f  =========\n', sh);
    fprintf('%-5s %11s %11s %11s %13s\n','','Original','PCA-based','EMD-based','DAGAF-based');
    rows = {'PE','SD','RNG'};
    for i = 1:3
        fprintf('%-5s%11.3f%11.3f%11.3f%13.3f\n', rows{i}, V4(i,1),V4(i,2),V4(i,3),V4(i,4));
        fprintf('%-5s%11.2f%11.2f%11.2f%13s   <- published\n','', PUB(i,1),PUB(i,2),PUB(i,3),'-');
    end
    % distance to the published row, first three columns only
    d = abs(V4(:,1:3) - PUB);
    fprintf('mean |diff| vs published:  PE %.3f   SD %.3f   RNG %.3f\n', ...
            mean(d(1,:)), mean(d(2,:)), mean(d(3,:)));
end

function img = drop(imdl,vh,vi,C,k)
w = ones(1,size(C,1)); w(k) = 0;
v = vi; v.meas = (w*C).';
img = inv_solve(imdl,vh,v); img.calc_colours.ref_level = 0;
end
