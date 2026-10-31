%% EIT_SHOW_PHANTOM  Render the FEM phantom used throughout this project.
% Geometry verbatim from "Error quantification_GREIT Test Parameters.m":
%   cylinder   ng_mk_cyl_models([2,2,0.2],[16,1],[0.1])  -- radius 2, height 2
%   lungs      sphere(+-1, 0.1, 1 ; r),  r = 0.3..0.6,  sigma = 0.5
%   heart      sphere(  0, 0.1, 1 ; 0.3),               sigma = 1.5
%   background sigma = 1
%
% show_fem draws unclipped lines when used inside subplot(), so each 3-D
% view is rendered to its own figure and the panels are composed afterwards
% with imread/imshow.

clc; close all;
SIG_L = 0.5; SIG_H = 2.0; RH = 0.3;
stim = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
tmp  = tempname; mkdir(tmp);

img03 = build(0.3, RH, SIG_L, SIG_H, stim);
img06 = build(0.6, RH, SIG_L, SIG_H, stim);

%% --- 3-D panels, one figure each ----------------------------------------
% show_fem draws stray slivers outside the axes that survive both clipping
% and rasterisation, so the 3-D views are drawn directly from the mesh with
% trisurf instead. All geometry is taken from the model itself.
% Panel (a) is drawn from the mesh built above by draw3d() and also written on
% its own as fig1a_phantom.png (used by eit_graphical_abstract). Until
% 2026-10-31 this panel imported a render first published as Fig. (2) of
% Cheng, Su and Ko 2022 (Curr Med Imaging 18:1396); that import and
% prep_fig1a.m were removed so that the repository holds only its own images.
p3d = fullfile(pwd,'fig1a_phantom.png');
draw3d(img06, p3d);

%% --- composite -----------------------------------------------------------
xax = -2:4/63:2;
f = figure('Color','w','Position',[0 0 1460 720],'Visible','off');

% --- top row: geometry ---------------------------------------------------
subplot(2,6,[1 2 3]); imshow(imread(p3d));
title(sprintf('(a) FEM model, 16 electrodes\ncylinder r=2, h=2, r_{lung}=0.6'),'FontSize',10);

subplot(2,6,[4 5 6]);
th = linspace(0,2*pi,400);
plot(2*cos(th), 2*sin(th), 'k-', 'LineWidth', 1.5); hold on; axis equal
fill(0.6*cos(th)+1, 0.6*sin(th)+0.1, [0.30 0.45 0.85], 'FaceAlpha',0.65, 'EdgeColor','none');
fill(0.6*cos(th)-1, 0.6*sin(th)+0.1, [0.30 0.45 0.85], 'FaceAlpha',0.65, 'EdgeColor','none');
fill(0.3*cos(th),   0.3*sin(th)+0.1, [0.85 0.25 0.25], 'FaceAlpha',0.75, 'EdgeColor','none');
E = img06.fwd_model.electrode;
for k = 1:numel(E)
    c = mean(img06.fwd_model.nodes(E(k).nodes,:),1);
    plot(c(1), c(2), 'o', 'MarkerSize',7, 'MarkerFaceColor',[0.1 0.6 0.2], ...
         'MarkerEdgeColor','k');
    text(c(1)*1.16, c(2)*1.16, num2str(k), 'FontSize',7.5, ...
         'HorizontalAlignment','center');
end
xlim([-2.5 2.5]); ylim([-2.5 2.5]); axis off
title({'(b) top view at the electrode plane','adjacent drive, 16\times13 = 208 measurements'},'FontSize',10);

% --- bottom row: conductivity -------------------------------------------
xax = -2:4/63:2;
slots = {[7 8],[9 10]};
for i = 1:2
    if i==1, im = img03; rl = 0.3; else, im = img06; rl = 0.6; end
    s = calc_slices(im);
    subplot(2,6,slots{i});
    imagesc(xax,xax,s,[0.4 2.1]); axis image; set(gca,'YDir','normal');
    colormap(gca,parula); colorbar; hold on
    plot(xax, 0.1*ones(size(xax)),'w:','LineWidth',1.2);
    title(sprintf('(%c) \\sigma at z=1,  r_{lung}=%.1f',char('b'+i),rl),'FontSize',10);
    xlabel('x'); ylabel('y');
end

subplot(2,6,[11 12]);
s3 = calc_slices(img03); s6 = calc_slices(img06);
plot(xax, s3(31,:),'LineWidth',1.8); hold on
plot(xax, s6(31,:),'LineWidth',1.8);
ylim([0.3 2.2]); grid on; axis square; xlabel('x'); ylabel('\sigma');
legend('r_{lung}=0.3','r_{lung}=0.6','Location','southeast'); legend boxoff
title('(e) horizontal centre cut (row 31)','FontSize',10);

exportgraphics(f,'fig_phantom.png','Resolution',200);
fprintf('wrote fig_phantom.png\n');
rmdir(tmp,'s');

function img = build(rl, RH, SIG_L, SIG_H, stim)
    lungs = sprintf('solid obj1 = sphere(1,0.1,1;%.2f) or sphere(-1,0.1,1;%.2f);',rl,rl);
    heart = sprintf('solid obj2 = sphere(0,0.1,1;%.2f);',RH);
    [f,midx] = ng_mk_cyl_models([2,2,0.2],[16,1],[0.1],{'obj1','obj2',[lungs heart]});
    f.stimulation = stim;
    img = mk_image(f,1);
    img.elem_data(midx{2}) = SIG_L;
    img.elem_data(midx{3}) = SIG_H;
end

function draw3d(img, outfile)
% Panel (a): the outer surface mesh as a wireframe, the lung and heart regions
% as solid surfaces, the electrodes as their own surface patches. Drawn from the
% model built here; no external render is used.
f = img.fwd_model;  N = f.nodes;
sig = img.elem_data;
isL = abs(sig-0.5) < 1e-9;  isH = abs(sig-2.0) < 1e-9;
h = figure('Color','w','Position',[0 0 900 760],'Visible','off');
ax = axes(h); hold(ax,'on');
trisurf(f.boundary, N(:,1),N(:,2),N(:,3), 'Parent',ax, ...
        'FaceColor','none','EdgeColor',[0.15 0.15 0.15],'EdgeAlpha',0.55,'LineWidth',0.35);
% inclusion surfaces: boundary faces of each region, kept only if every node lies
% on that sphere (find_boundary also returns a few faces spanning the region)
sph = {[ 1 0.1 1 0.6], [-1 0.1 1 0.6], [0 0.1 1 0.3]};
col = {[0.12 0.40 0.95], [0.12 0.40 0.95], [0.60 0.06 0.06]};
reg = {isL, isL, isH};
for s = 1:3
    F  = find_boundary(f.elems(reg{s},:));
    c  = sph{s}(1:3);  r = sph{s}(4);
    d  = sqrt(sum((N - c).^2, 2));
    F  = F(all(abs(d(F) - r) < 0.12*r + 0.02, 2), :);
    % faces only, then the facet edges as explicit line segments: the software
    % renderer used in batch mode draws a few heart edges of a patch as long
    % stray rays, explicit segments do not
    patch(ax,'Faces',F,'Vertices',N,'FaceColor',col{s},'EdgeColor','none', ...
          'FaceLighting','flat','AmbientStrength',0.55);
    E = unique(sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])],2),'rows');
    X = [N(E(:,1),1) N(E(:,2),1) nan(size(E,1),1)].';
    Y = [N(E(:,1),2) N(E(:,2),2) nan(size(E,1),1)].';
    Z = [N(E(:,1),3) N(E(:,2),3) nan(size(E,1),1)].';
    % lift the edges a hair off the surface toward the sphere's outside so they are not hidden
    plot3(ax, c(1)+(X-c(1))*1.003, c(2)+(Y-c(2))*1.003, c(3)+(Z-c(3))*1.003, '-', ...
          'Color',col{s}*0.55,'LineWidth',0.4);
end
for k = 1:numel(f.electrode)
    en = f.electrode(k).nodes;
    bf = f.boundary(all(ismember(f.boundary, en),2),:);
    if isempty(bf)
        c = mean(N(en,:),1);
        plot3(ax,c(1),c(2),c(3),'o','MarkerSize',8,'MarkerFaceColor',[0.05 0.45 0.10],'MarkerEdgeColor','k');
    else
        trisurf(bf, N(:,1),N(:,2),N(:,3),'Parent',ax,'FaceColor',[0.05 0.45 0.10],'EdgeColor','none');
    end
end
axis(ax,'equal'); view(ax,-37.5,30);
xlim(ax,[-2 2]); ylim(ax,[-2 2]); zlim(ax,[0 2]);
grid(ax,'off'); box(ax,'off');
set(ax,'FontSize',11,'XTick',-2:0.5:2,'YTick',-2:0.5:2,'ZTick',0:0.5:2);
camlight(ax,'headlight'); material(ax,'dull');
exportgraphics(h,outfile,'Resolution',200); close(h);
fprintf('wrote %s\n', outfile);
end
