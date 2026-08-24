%% EIT_AXIS_TEST  Does filtering along the injection axis reduce ringing?
% Reference points (heart target, r_lung 0.6, sigma_heart 1.5):
%   Original       SD 0.981  RNG 38.409
%   PCA  #4        SD 0.025  RNG  1.464
%   EMD  #1        SD 0.059  RNG  1.088
%   DAGAF chan cap6 #1  SD 0.030  RNG  7.768
%   DAGAF chan cap4 #1  SD 0.046  RNG  4.939
clc;
stim = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
lungs='solid obj1 = sphere(1,0.1,1;0.60) or sphere(-1,0.1,1;0.60);';
heart='solid obj2 = sphere(0,0.1,1;0.30);';
[fm,midx]=ng_mk_cyl_models([2,2,0.2],[16,1],[0.1],{'obj1','obj2',[lungs heart]});
fm.stimulation=stim; im=mk_image(fm,1); vh=fwd_solve(im);
im.elem_data(midx{2})=0.5; im.elem_data(midx{3})=1.5; vi=fwd_solve(im);
X=double(vi.meas(:)).';
opt=struct('imgsz',[64 64],'distr',0,'Nsim',1000,'target_size',0.02,'noise_figure',0.5);
fg=ng_mk_cyl_models([2,2,0.2],[16,1],[0.1]); fg.stimulation=stim; fg=mdl_normalize(fg,1);
imdl=mk_GREIT_model(mk_image(fg,1),0.2,[],opt);
tgt=[0;0.1;1;0.30];

fprintf('\n%-6s %-5s %-6s %6s %8s %8s %9s\n','axis','cap','drop','M','PE','SD','RNG');
best=struct('rng',inf);
for cap = [1 2 3 4 5 7]
    [C,info] = dagaf_axis(X,'inj',2,2,'p',40,1e-8,3,cap);
    err = norm(sum(C,1)-X)/norm(X);
    if err > 1e-10, fprintf('  completeness FAILED (%.1e)\n', err); continue; end
    Mmed = median(info.M(:,1));
    for k = 1:size(C,1)-1
        w=ones(1,size(C,1)); w(k)=0;
        v=vi; v.meas=(w*C).';
        p=eval_GREIT_fig_merit(inv_solve(imdl,vh,v),tgt);
        star='';
        if p(4)<0.10 && p(5)<best.rng
            best=struct('rng',p(5),'sd',p(4),'pe',p(2),'cap',cap,'k',k,'M',Mmed);
            star='  *';
        end
        fprintf('%-6s %-5d %-6d %6d %8.3f %8.3f %9.3f%s\n','inj',cap,k,Mmed,p(2),p(4),p(5),star);
    end
end
if isfinite(best.rng)
    fprintf('\nbest (SD<0.10): cap=%d drop #%d  M=%d -> PE %.3f SD %.3f RNG %.3f\n', ...
            best.cap,best.k,best.M,best.pe,best.sd,best.rng);
    fprintf('vs DAGAF chan cap4: SD 0.046 RNG 4.939   |   PCA: SD 0.025 RNG 1.464\n');
else
    fprintf('\nno configuration reached SD < 0.10 on the injection axis\n');
end
