%% EIT_DAGAF_RINGING
% Attack the ringing (RNG) weakness of DAGAF-based heart imaging.
% Hypothesis: IMF1 with a 13-point window removes the sharp peaks at the
% injection-block boundaries wholesale, taking the heart's high spatial
% frequency content with it. A narrower window, or a different boundary
% extension, should trade SD against RNG differently.
%
% Sweeps maskCap x eType x chi, always dropping IMF1, heart target.
% Reference points: Original SD 0.981 RNG 38.409 | PCA#4 SD 0.025 RNG 1.464
% Run setup_eit_env first.

clc;
P.r_lung = 0.6; P.r_heart = 0.3;
stim = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
lungs = sprintf('solid obj1 = sphere(1,0.1,1;%.2f) or sphere(-1,0.1,1;%.2f);',P.r_lung,P.r_lung);
heart = sprintf('solid obj2 = sphere(0,0.1,1;%.2f);',P.r_heart);
[fmdl,midx] = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1],{'obj1','obj2',[lungs heart]});
fmdl.stimulation = stim;
im = mk_image(fmdl,1); vh = fwd_solve(im);
im.elem_data(midx{2}) = 0.5; im.elem_data(midx{3}) = 2.0; vi = fwd_solve(im);
X = double(vi.meas(:)).';

opt = struct('imgsz',[64 64],'distr',0,'Nsim',1000,'target_size',0.02,'noise_figure',0.5);
fg = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1]); fg.stimulation = stim; fg = mdl_normalize(fg,1);
imdl = mk_GREIT_model(mk_image(fg,1),0.2,[],opt);
tgt = [0;0.1;1;P.r_heart];

p0 = eval_GREIT_fig_merit(inv_solve(imdl,vh,vi), tgt);
fprintf('\nreference   Original   PE %6.3f  SD %6.3f  RNG %8.3f\n', p0(2),p0(4),p0(5));
fprintf('reference   PCA #4     PE  0.049  SD  0.025  RNG    1.464\n\n');

caps   = [2 3 4 6 8 10 13 Inf];
etypes = {'p','d','r','c'};
chis   = [1.2 2 3];

fprintf('%7s %6s %5s %6s %8s %8s %9s\n','maskCap','eType','chi','M','PE','SD','RNG');
best = struct('rng',inf);
res = [];
for ci = 1:numel(chis)
  for ei = 1:numel(etypes)
    for k = 1:numel(caps)
        cap = caps(k);
        [di,dr,info] = dagaf2(X, 2, chis(ci), etypes{ei}, etypes{ei}, ...
                              40, 1e-8, 3, cap);
        if isempty(di), continue; end
        K = size(di,1); w = ones(1,K+1); w(1) = 0;
        v = vi; v.meas = (w*[di;dr]).';
        p = eval_GREIT_fig_merit(inv_solve(imdl,vh,v), tgt);
        M1 = info.M(1);
        capstr = 'Inf'; if isfinite(cap), capstr = num2str(cap); end
        star = '';
        % only count it as an improvement if the shape stays good
        if p(4) < 0.10 && p(5) < best.rng
            best = struct('rng',p(5),'sd',p(4),'pe',p(2),'cap',capstr, ...
                          'e',etypes{ei},'chi',chis(ci),'M',M1);
            star = '  *';
        end
        fprintf('%7s %6s %5.1f %6d %8.3f %8.3f %9.3f%s\n', ...
                capstr, etypes{ei}, chis(ci), M1, p(2), p(4), p(5), star);
        res(end+1,:) = [cap chis(ci) ei M1 p(2) p(4) p(5)]; %#ok<SAGROW>
    end
  end
end
fprintf('\nbest (SD < 0.10, lowest RNG): maskCap=%s eType=%s chi=%.1f M=%d\n', ...
        best.cap, best.e, best.chi, best.M);
fprintf('   -> PE %.3f  SD %.3f  RNG %.3f\n', best.pe, best.sd, best.rng);
save('dagaf_ringing_sweep.mat','res','caps','etypes','chis');
