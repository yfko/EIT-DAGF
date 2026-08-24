%% PREP_FIG1A  Prepare the panel (a) phantom render supplied by the authors.
% Reads the 3-D view exported from the Graphic Table of Contents deck, drops
% the embedded title band, trims the white margin, and writes the PNG that
% eit_show_phantom.m imports as panel (a).
A = imread('fig1a_source.tif');
if size(A,3) == 1, A = repmat(A,1,1,3); end
if ~isa(A,'uint8')
    A = uint8(double(A) / double(max(A(:))) * 255);   % no Image Processing Toolbox here
end
fprintf('source %d x %d\n', size(A,1), size(A,2));

% drop the title band: the plot content begins well below it
A = A(round(0.12*size(A,1)):end, :, :);

% trim the white margin
g   = min(A,[],3);
rows = find(any(g < 245, 2));
cols = find(any(g < 245, 1));
pad  = 12;
r1 = max(1,rows(1)-pad);  r2 = min(size(A,1), rows(end)+pad);
c1 = max(1,cols(1)-pad);  c2 = min(size(A,2), cols(end)+pad);
A = A(r1:r2, c1:c2, :);
fprintf('cropped %d x %d\n', size(A,1), size(A,2));
imwrite(A, 'fig1a_phantom.png');
