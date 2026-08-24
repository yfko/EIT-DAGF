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
p3d = fullfile(tmp,'a.png');
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

exportgraphics(f,'fig_phantom.png','Resolution',150);
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
f = img.fwd_model;  N = f.nodes;
sig = img.elem_data;
isL  = abs(sig-0.5) < 1e-9;                % lungs
isH  = abs(sig-2.0) < 1e-9;                % heart

h = figure('Color','w','Position',[0 0 760 660],'Visible','off');
% outer surface, translucent
% Use the boundary Netgen/EIDORS stored: find_boundary() returns fewer faces
% (5394 vs 5988) and draws stray edges across the axes.
trisurf(f.boundary, N(:,1),N(:,2),N(:,3), ...
        'FaceColor',[0.82 0.84 0.88],'FaceAlpha',0.16, ...
        'EdgeColor',[0.45 0.45 0.5],'EdgeAlpha',0.20); hold on
% lungs and heart
trisurf(find_boundary(f.elems(isL,:)), N(:,1),N(:,2),N(:,3), ...
        'FaceColor',[0.20 0.32 0.80],'FaceAlpha',0.95,'EdgeColor','none');
trisurf(find_boundary(f.elems(isH,:)), N(:,1),N(:,2),N(:,3), ...
        'FaceColor',[0.85 0.18 0.18],'FaceAlpha',0.98,'EdgeColor','none');
% electrodes
for k = 1:numel(f.electrode)
    c = mean(N(f.electrode(k).nodes,:),1);
    plot3(c(1),c(2),c(3),'o','MarkerSize',8, ...
          'MarkerFaceColor',[0.10 0.62 0.20],'MarkerEdgeColor','k');
end
axis equal off; view(-37,20); camlight headlight; lighting gouraud; material dull
xlim([-2.1 2.1]); ylim([-2.1 2.1]); zlim([-0.05 2.05]);
exportgraphics(h,outfile,'Resolution',200,'ContentType','image'); close(h);
end
