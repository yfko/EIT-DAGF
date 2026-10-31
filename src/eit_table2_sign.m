function eit_table2_sign()
%EIT_TABLE2_SIGN  Signed Pearson correlations behind Table 2 (round 10, panel item R17).
% Table 2 of the manuscript reports |r| of each component with the lung-only and
% heart-only measurement signatures (eit_revision_r1.m stores abs values only).
% This script recomputes the same correlations with their sign, on the same frame
% and signatures (the companion study's cached environment p2_env.mat, read only;
% its vendor/ copy of this code is at tag v1.0). No reconstruction is involved.
% Deposit copy: the companion path is resolved by companion_root (end of file).
% Writes table2_sign.txt (diary) next to this file.
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'rcada_emd'));
out = fullfile(here,'table2_sign.txt'); if isfile(out), delete(out); end
diary(out);
fprintf('eit_table2_sign  %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
S = load(fullfile(companion_root(here),'results','p2_env.mat'),'E');
E = S.E;
X  = reshape(double(E.X),1,[]);
dl = reshape(double(E.d_lung),1,[]);
dh = reshape(double(E.d_heart),1,[]);
fprintf('||d_lung|| %.5g  ||d_heart|| %.5g  corr %.3f  (revision_r1_results.txt: 0.02457, 0.0007319, 0.253)\n', ...
        norm(dl), norm(dh), pear(dl,dh));
% --- PCA as in eit_revision_r1.m / p2_decompose.m
A = reshape(X(:),13,16); mu = mean(A,1); Ac = A - repmat(mu,13,1);
[V,Dg] = eig(cov(Ac)); [~,ord] = sort(diag(Dg),'descend'); V = V(:,ord); PC = Ac*V;
cp = zeros(17,208); for k = 1:16, cp(k,:) = reshape(PC(:,k)*V(:,k)',1,[]); end
cp(17,:) = reshape(repmat(mu,13,1),1,[]);
% --- EMD (RCADA eemd, Nstd 0, NE 1) and DAGAF (cap 6) as in eit_compare_table1.m
[imf,res] = emd_rcada(X);  ce = [imf; res];
[imf,res] = dagaf2(X, 6, 2, 'p', 'p', 40, 1e-8, 3, 6);  cd_ = [imf; res];
D = {'PCA', cp, 4; 'EMD', ce, 1; 'DAGAF', cd_, 1};
for m = 1:3
    C = D{m,2}; nk = size(C,1) - 1;
    fprintf('\n%s (closure %.2g)\n', D{m,1}, norm(sum(C,1) - X)/norm(X));
    fprintf('  %-28s', 'signed r with lung signature'); fprintf(' %8.3f', arrayfun(@(k) pear(C(k,:),dl), 1:nk)); fprintf('\n');
    fprintf('  %-28s', 'signed r with heart signature'); fprintf(' %8.3f', arrayfun(@(k) pear(C(k,:),dh), 1:nk)); fprintf('\n');
    k = D{m,3};
    fprintf('  separating component #%d: r_lung %+.3f, r_heart %+.3f;  most lung-correlated #2: r_lung %+.3f, r_heart %+.3f\n', ...
            k, pear(C(k,:),dl), pear(C(k,:),dh), pear(C(2,:),dl), pear(C(2,:),dh));
    % raw (uncentred) projection of the separating component on the heart signature: C_k = a*d_heart + rest
    a = (C(k,:)*dh.')/(dh*dh.');
    fprintf('  projection coefficient a = <C_%d, d_heart>/<d_heart, d_heart> = %+.4g;  ||a d_heart||/||C_%d|| = %.3f\n', k, a, k, norm(a*dh)/norm(C(k,:)));
end
eit_table2_sign_images();
diary off
end
function r = pear(a,b)
R = corrcoef(a(:), b(:)); r = R(1,2);
end

function eit_table2_sign_images()
%EIT_TABLE2_SIGN_IMAGES  Image-domain companion to eit_table2_sign (R17): under the
% reference GREIT model of the companion environment, reconstruct as difference data the
% heart-only signature d_heart and each separating component C_j, and report the heart-region
% sign of recon(d_heart), of -recon(C_j) (the removed-component term) and their correlation.
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'rcada_emd'));
S = load(fullfile(companion_root(here),'results','p2_env.mat'),'E');  E = S.E;
X  = reshape(double(E.X),1,[]);  vh = double(E.vh.meas(:)).';
roi = logical(E.roi_h(:));
rec = @(d) local_rec(E, vh, d);
eh = rec(reshape(double(E.d_heart),1,[]));
fprintf('\nrecon(d_heart): heart-region mean %+.4g, image mean %+.4g, heart-region share of |image| %.2f\n', ...
        mean(eh(roi)), mean(eh), sum(abs(eh(roi)))/sum(abs(eh)));
A = reshape(X(:),13,16); mu = mean(A,1); Ac = A - repmat(mu,13,1);
[V,Dg] = eig(cov(Ac)); [~,ord] = sort(diag(Dg),'descend'); V = V(:,ord); PC = Ac*V;
c4 = reshape(PC(:,4)*V(:,4)',1,[]);
[imf,~] = emd_rcada(X);  c_emd = imf(1,:);
[imf,~] = dagaf2(X, 6, 2, 'p', 'p', 40, 1e-8, 3, 6);  c_dag = imf(1,:);
D = {'PCA PC4', c4; 'EMD mode 1', c_emd; 'DAGAF mode 1', c_dag};
for m = 1:3
    t = -rec(D{m,2});          % removed-component term
    fprintf('%-13s removed-component term -recon(C_j): heart-region mean %+.4g; r_h(term, recon(d_heart)) %+.3f; r_h(term, -recon(d_heart)) %+.3f\n', ...
            D{m,1}, mean(t(roi)), pear(t(roi), eh(roi)), pear(t(roi), -eh(roi)));
end
end
function e = local_rec(E, vh, d)
v = E.vi;  v.meas = (vh + d).';
img = inv_solve(E.imdl, E.vh, v);  e = img.elem_data;
end

function p = companion_root(here)
%COMPANION_ROOT  Locate the companion repository EIT_noise_geometry_pre-registered.
% Set the environment variable EIT_NOISE_GEOMETRY to its root, or clone it next
% to this repository (../EIT_noise_geometry_pre-registered relative to the
% repository root). https://github.com/yfko/EIT_noise_geometry_pre-registered, tag v1.0.
p = getenv('EIT_NOISE_GEOMETRY');
if isempty(p), p = fullfile(here,'..','..','EIT_noise_geometry_pre-registered'); end
if ~isfolder(fullfile(p,'results'))
    error('companion_root:notFound', ['Companion repository not found at %s. Clone ' ...
        'https://github.com/yfko/EIT_noise_geometry_pre-registered (tag v1.0) or set EIT_NOISE_GEOMETRY.'], p);
end
end
