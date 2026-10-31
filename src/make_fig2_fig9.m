function make_fig2_fig9()
%MAKE_FIG2_FIG9  Assemble manuscript Figures 2 and 9 from existing results.
%
% Deposit copy (github.com/yfko/EIT-DAGF). Read-only on the companion repository
% EIT_noise_geometry_pre-registered (tag v1.0): loads its results/*.mat
% (including the cached reference GREIT model results/p2_env.mat) and writes
% nothing there. Writes only into the folder that holds this script:
%   fig2_split_composite.png   (300 dpi)
%   fig9_window_sweep.png      (300 dpi)
%   make_fig2_fig9_log.txt     (consistency checks)
%
% Display convention as the companion's p2_split_figure.m: calc_slices row 1 is
% max y, imagesc's default reverse YDir puts it at the top -> anterior up.
% Colour map: symmetric blue-white-red, limits [-c c].

here = fileparts(mfilename('fullpath'));
P2   = companion_root(here);                 % deposit copy: see companion_root at the end of this file
R    = fullfile(P2,'results');
P1   = here;                                 % human_separating_component.mat, written here by eit_human_separating_component
logf = fullfile(here,'make_fig2_fig9_log.txt');
if isfile(logf), delete(logf); end
diary(logf);
fprintf('make_fig2_fig9  %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>

if exist('calc_slices','file') ~= 2
    error('make_fig2_fig9:eidors', 'Start EIDORS first (setup_eit_env or run its startup.m).');
end

%% ===================== data: phantom, reference position ====================
A = load(fullfile(R,'p2_separating_component.mat'));      SC = A.S;
B = load(fullfile(R,'p2_separating_component_ref.mat'));  T  = B.T;
C = load(fullfile(R,'p2_separating_component_offctr.mat')); O = C.S;
V = load(fullfile(R,'p2_env.mat'));                        E  = V.E;

% image struct on the reference GREIT model, for slicing element data
img0 = inv_solve(E.imdl, E.vh, E.vi);  img0.calc_colours.ref_level = 0;
chk = norm(img0.elem_data - SC.unproc.elem) / norm(SC.unproc.elem);
fprintf('reference GREIT model reproduces the saved unprocessed image: rel diff %.2e\n', chk);
assert(chk < 1e-10, 'p2_env.mat model does not match p2_separating_component.mat');
sl = @(img, e) squeeze(calc_slices(setfield(img,'elem_data',e))); %#ok<SFLD>

refNames = {'PCA','EMD','DAGAF','HP-gauss-chan'};
refLab   = {'PCA (PC4)','EMD (mode 1)','DAGAF (mode 1)','Gaussian control'};
rows = struct('grp',{},'lab',{},'p',{},'o',{},'a',{},'r',{});
for i = 1:numel(refNames)
    m = find(strcmp(SC.names, refNames{i}));
    t = find(strcmp({T.m.name}, refNames{i}));
    assert(isscalar(m) && isscalar(t));
    rows(end+1) = struct('grp','Phantom, reference position (0, 0.1)', 'lab',refLab{i}, ... %#ok<AGROW>
        'p',SC.slices.unproc, 'o',SC.slices.out{m}, 'a',SC.slices.asep{m}, ...
        'r',sl(img0, -T.m(t).e_sep_ref));
    fprintf('ref  %-14s key %d   max|out| %.4g  max|-sep| %.4g  max|ref part| %.4g  max|phys| %.4g\n', ...
        refNames{i}, SC.m(m).key, max(abs(rows(end).o(:))), max(abs(rows(end).a(:))), ...
        max(abs(rows(end).r(:))), max(abs(rows(end).p(:))));
end

%% ===================== data: phantom, heart at (0, -0.401) ==================
ip = find(arrayfun(@(p) isequal(round(p.heart_ctr,3), [0 -0.401]), O.pos));
assert(isscalar(ip), 'position (0,-0.401) not found');
Q = O.pos(ip).Q;  ph = O.pos(ip).phys;
offNames = {'PCA','EMD','DAGAF'};  offLab = {'PCA (PC4)','EMD (mode 1)','DAGAF (mode 1)'};
for i = 1:3
    m = find(strcmp({Q.name}, offNames{i}));
    rows(end+1) = struct('grp','Phantom, heart at (0, -0.401)', 'lab',offLab{i}, ... %#ok<AGROW>
        'p',ph.slice, 'o',Q(m).slice_out, 'a',Q(m).slice_asep, 'r',[]);
    fprintf('off  %-14s key %d   max|out| %.4g  max|-sep| %.4g  max|phys| %.4g  split closure %.1e\n', ...
        offNames{i}, Q(m).key, max(abs(Q(m).slice_out(:))), max(abs(Q(m).slice_asep(:))), ...
        max(abs(ph.slice(:))), Q(m).closure_split);
end

%% ===================== data: human breath-hold frame 10 =====================
H = load(fullfile(P1,'human_separating_component.mat'));
k = find(H.frames == 10);  assert(k == 9);
% GREIT model exactly as eit_human_separating_component.m lines 32-43
fmdl = mk_library_model('adult_male_16el_lungs');
fmdl.electrode = fmdl.electrode([9:16,1:8]);
img = mk_image(fmdl, 1);
img.elem_data([fmdl.mat_idx{2}; fmdl.mat_idx{3}]) = 0.3;
[stim,msel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
img.fwd_model.stimulation = stim;
img.fwd_model = mdl_normalize(img.fwd_model, 1);
opt = struct('imgsz',[32 32],'distr',3,'Nsim',500,'target_size',0.03, ...
             'target_offset',0,'noise_figure',0.5,'square_pixels',1);
imdlH = mk_GREIT_model(img, 0.25, [], opt);
imdlH.fwd_model.meas_select = msel;
D = load('montreal_data_1995');  Z = double(D.zc_breathhold);
imH = inv_solve(imdlH, Z(:,1), Z(:,10));  imH.calc_colours.ref_level = 0;
chk = norm(imH.elem_data - H.PCA.e_phys(:,k)) / norm(H.PCA.e_phys(:,k));
fprintf('rebuilt human GREIT model reproduces the saved frame-10 physiological term: rel diff %.2e\n', chk);
assert(chk < 1e-6, 'rebuilt human model does not match human_separating_component.mat');
humLab = {'PCA (PC3)','EMD (mode 1)','DAGAF (mode 1)'};
for i = 1:3
    Rm = H.(H.methods{i});
    rows(end+1) = struct('grp','Human breath-hold, frame 10', 'lab',humLab{i}, ... %#ok<AGROW>
        'p',sl(imH, Rm.e_phys(:,k)), 'o',sl(imH, Rm.e_out(:,k)), 'a',sl(imH, Rm.e_art(:,k)), 'r',[]);
    fprintf('hum  %-14s max|out| %.4g  max|-rem| %.4g  max|phys| %.4g  |out-(phys+art)|/|out| %.1e\n', ...
        H.methods{i}, max(abs(rows(end).o(:))), max(abs(rows(end).a(:))), max(abs(rows(end).p(:))), ...
        norm(Rm.e_out(:,k)-(Rm.e_phys(:,k)+Rm.e_art(:,k)))/norm(Rm.e_out(:,k)));
end

%% ===================== Figure 2 =============================================
cols = {'physiological term','contract output', ...
        [char(8722) 'recon(removed)'],'reference-frame part'};
nr = numel(rows);  nc = 4;
cmap = bluewhitered(256);
W = 7.0;  ph_ = 1.02;  top = 0.42;  gapG = 0.22;  bottom = 0.38;
grpBreak = [0, find(~strcmp({rows(2:end).grp}, {rows(1:end-1).grp}))];  % index before each new group
nG = numel(grpBreak);
Hh = top + nr*ph_ + (nG-1)*gapG + bottom;
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 W Hh],'Visible','off');
colormap(fig, cmap);
x0 = 1.30;  pw = 1.30;  dx = 1.40;
y = Hh - top;
for r = 1:nr
    if r > 1 && ~strcmp(rows(r).grp, rows(r-1).grp), y = y - gapG; end
    if r == 1 || ~strcmp(rows(r).grp, rows(r-1).grp)
        annotation(fig,'textbox',[0.005, (y+0.03)/Hh, 0.99, 0.18/Hh],'Units','normalized', ...
            'String',rows(r).grp,'EdgeColor','none','FontSize',8.5,'FontWeight','bold', ...
            'VerticalAlignment','bottom','Interpreter','none');
    end
    y = y - ph_;
    cm = max(abs(rows(r).o(:)));
    panels = {rows(r).p, rows(r).o, rows(r).a, rows(r).r};
    for c = 1:nc
        ax = axes(fig,'Units','inches','Position',[x0+(c-1)*dx, y+0.04, pw, ph_-0.08]);
        M = panels{c};
        if isempty(M)
            axis(ax,'off');   % caption: blank for rows other than the phantom reference position
            continue
        end
        h = imagesc(ax, M);  set(h,'AlphaData',~isnan(M));
        axis(ax,'image','off');  clim(ax,[-cm cm]);
        if r == 1, title(ax, cols{c}, 'FontSize',8.5,'FontWeight','normal'); end
    end
    annotation(fig,'textbox',[0.005, y/Hh, (x0-0.05)/W, ph_/Hh], 'String', ...
        {rows(r).lab, sprintf('scale \\pm%.0f', cm)}, 'EdgeColor','none', ...
        'FontSize',7.5,'HorizontalAlignment','right','VerticalAlignment','middle');
end
% colour bar (shared shape; limits are per row)
cb = axes(fig,'Units','inches','Position',[x0+0.2, 0.20, 3*dx-0.4, 0.07]);
imagesc(cb, linspace(-1,1,256)); colormap(cb, cmap); set(cb,'YTick',[],'XTick',[1 128.5 256], ...
    'XTickLabel',{[char(8722) 'scale'],'0','+scale'},'FontSize',7,'TickLength',[0 0]);
f2 = fullfile(here,'fig2_split_composite.png');
exportgraphics(fig, f2, 'Resolution',300);  close(fig);
fprintf('wrote %s\n', f2);

%% ===================== Figure 9 =============================================
Wn = load(fullfile(R,'p2_separating_component_window.mat'));  S9 = Wn.S;
L = S9.L;  kinds = {'gauss','box'};  kname = {'Gaussian','rectangular'};
Lshow = [7 9 13 17];
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 7.2 3.3],'Visible','off');
ax = axes(fig,'Units','inches','Position',[0.55 0.50 2.35 2.45]);
hold(ax,'on');  mk = {'-o','-s'};  colr = [0.15 0.35 0.70; 0.80 0.35 0.10];
for ik = 1:2
    sel = strcmp({S9.r.kind}, kinds{ik});
    LL = [S9.r(sel).L];  cc = [S9.r(sel).conc];
    plot(ax, LL, cc, mk{ik}, 'Color',colr(ik,:),'MarkerFaceColor',colr(ik,:), ...
         'MarkerSize',4,'LineWidth',1.2,'DisplayName',kname{ik});
    fprintf('conc %-6s: %s\n', kinds{ik}, sprintf('L=%d %.2f  ', [LL; cc]));
end
xline(ax, 13, '--', 'block length 13', 'Color',[.45 .45 .45], 'LineWidth',0.9, 'FontSize',7, ...
      'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','left', 'HandleVisibility','off');
set(ax,'XTick',L,'FontSize',7.5,'Box','off');  xlim(ax,[2 30]);  ylim(ax,[0 5.5]);
xlabel(ax,'window length L (samples)','FontSize',8);
ylabel(ax,'concentration (RMS heart / RMS image)','FontSize',8);
legend(ax,'Location','north','NumColumns',2,'Box','off','FontSize',7);
title(ax,'(a)','FontSize',8.5,'FontWeight','normal');
colormap(fig, cmap);
for ik = 1:2
    for il = 1:numel(Lshow)
        j = find(strcmp({S9.r.kind}, kinds{ik}) & [S9.r.L] == Lshow(il));
        rr = S9.r(j);  M = rr.slice;  cm = max(abs(M(:)));
        a2 = axes(fig,'Units','inches','Position',[3.25+(il-1)*0.98, 1.72-(ik-1)*1.45, 0.92, 0.92]);
        h = imagesc(a2, M);  set(h,'AlphaData',~isnan(M));
        axis(a2,'image','off');  clim(a2,[-cm cm]);
        title(a2, sprintf('L = %d, conc %.2f', rr.L, rr.conc), 'FontSize',7,'FontWeight','normal');
        text(a2, 0.5, -0.06, sprintf('\\pm%.0f', cm), 'Units','normalized', ...
             'HorizontalAlignment','center','FontSize',6.5,'Color',[.4 .4 .4]);
        if il == 1
            text(a2, -0.12, 0.5, kname{ik}, 'Units','normalized','Rotation',90, ...
                 'HorizontalAlignment','center','FontSize',8);
        end
    end
end
annotation(fig,'textbox',[3.2/7.2, 3.02/3.3, 4.0/7.2, 0.25/3.3],'String', ...
    ['(b)  ' char(8722) 'recon(hp(v_h)), each on its own scale'],'EdgeColor','none', ...
    'HorizontalAlignment','center','FontSize',8.5);
for j = 1:numel(S9.r)
    M = S9.r(j).slice;  ctr = M(31:34,31:34);   % 4x4 block at the image centre (x=0, y~0)
    fprintf('window %-5s L=%2d  conc %.2f  centre mean / max|image| %+.3f\n', S9.r(j).kind, S9.r(j).L, ...
        S9.r(j).conc, mean(ctr(:))/max(abs(M(:))));
end
f9 = fullfile(here,'fig9_window_sweep.png');
exportgraphics(fig, f9, 'Resolution',300);  close(fig);
fprintf('wrote %s\n', f9);
diary off
end

% -------------------------------------------------------------------------
function s = num3(x)
s = sprintf('%.3g', x);
end
function c = bluewhitered(n)
h = n/2;  t = linspace(0,1,h).';
c = [[t, t, ones(h,1)]; [ones(h,1), flipud(t), flipud(t)]];
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
