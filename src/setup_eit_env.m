function setup_eit_env(eidors_dir)
%SETUP_EIT_ENV  Point MATLAB at EIDORS and the Netgen binary on this Mac.
%
%   setup_eit_env                       % looks for EIDORS in the usual places
%   setup_eit_env('/path/to/eidors')    % explicit
%
% Netgen was installed to /Applications/Netgen.app. EIDORS shells out to the
% `netgen` executable inside it, so both NETGENDIR and PATH must be set
% BEFORE eidors startup runs.
%
% NOTE: that Netgen build (6.1, 2016) is x86_64 only. On Apple Silicon it
% needs Rosetta 2:   softwareupdate --install-rosetta --agree-to-license

NG = '/Applications/Netgen.app/Contents/MacOS';

%% --- Netgen ---------------------------------------------------------------
if ~exist(fullfile(NG,'netgen'), 'file')
    error('netgen not found at %s -- is Netgen.app installed?', NG);
end
setenv('NETGENDIR', NG);
setenv('PATH', [NG ':' getenv('PATH')]);

[st, out] = system('netgen -batchmode -help 2>&1 | head -2');
if contains(out, 'Bad CPU type')
    error(['Netgen cannot run: Rosetta 2 is not installed.\n' ...
           'Run this in Terminal, then retry:\n' ...
           '    softwareupdate --install-rosetta --agree-to-license']);
end
fprintf('Netgen  : %s  (exit %d)\n', NG, st);

%% --- EIDORS ---------------------------------------------------------------
if nargin < 1 || isempty(eidors_dir)
    cand = {fullfile(getenv('HOME'),'Documents','MATLAB','eidors-v3.8','eidors'), ...
            fullfile(getenv('HOME'),'Documents','MATLAB','eidors','eidors'), ...
            fullfile(getenv('HOME'),'Documents','MATLAB','eidors'), ...
            fullfile(getenv('HOME'),'eidors'), ...
            fullfile(pwd,'eidors')};
    eidors_dir = '';
    for i = 1:numel(cand)
        if exist(fullfile(cand{i},'startup.m'), 'file')
            eidors_dir = cand{i}; break
        end
    end
    if isempty(eidors_dir)
        error(['EIDORS not found. Download eidors-v3.x from ' ...
               'https://sourceforge.net/projects/eidors3d/ , unzip it, ' ...
               'then call setup_eit_env(''/path/to/eidors'').']);
    end
end
run(fullfile(eidors_dir, 'startup.m'));
fprintf('EIDORS  : %s\n', eidors_dir);

%% --- this project ---------------------------------------------------------
here = fileparts(mfilename('fullpath'));
addpath(here);
fprintf('project : %s\n', here);

%% --- smoke test -----------------------------------------------------------
fprintf('\nsmoke test: building a 16-electrode cylinder via Netgen...\n');
fmdl = ng_mk_cyl_models([2,2,0.2], [16,1], [0.1]);
fprintf('  OK: %d nodes, %d elements\n', ...
        size(fmdl.nodes,1), size(fmdl.elems,1));

fprintf('\nsmoke test: dagaf2 on a 208-point vector...\n');
x = sin(2*pi*(0:207)/13) + 0.3*sin(2*pi*(0:207)/60);
[imf, res] = dagaf2(x, 6, 2, 'p', 'p', 40, 1e-8, 3, 6);
fprintf('  OK: %d IMFs, completeness err = %.2e\n', ...
        size(imf,1), norm(sum(imf,1)+res - x)/norm(x));

fprintf('\nReady. Now run:  eit_compare_table1\n');
end
