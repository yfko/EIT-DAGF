function [imf, res] = emd_rcada(x, ~)
%EMD_RCADA  Row-oriented EMD wrapper around the RCADA/NCU eemd.m package.
%
% The package (Zhaohua Wu, Research Center for Adaptive Data Analysis,
% National Central University) was distributed from rcada.ncu.edu.tw, which
% is now offline; these files were recovered from the Internet Archive
% snapshot of 2009-06-21 and live in ./rcada_emd.
%
%   allmode = eemd(Y, Nstd, NE)     Nstd=0, NE=1 degenerates to plain EMD.
%   allmode is N x (TNM+2) with TNM = fix(log2(N))-1:
%       column 1        = the original signal
%       columns 2..TNM+1 = IMFs, high to low frequency
%       column TNM+2    = residual
%   For N = 208, TNM = 6 -- which is exactly the "six IMFs" reported in
%   Cheng, Su & Ko 2022.
%
% This wrapper returns the IMFs as ROWS plus the residual separately, so it
% matches the interface used by dagaf2 and by the CMI paper's own scripts.

x = reshape(double(x), 1, []);
A = eemd(x(:), 0, 1);                 % N x (TNM+2)
imf = A(:, 2:end-1).';                % IMFs as rows
res = A(:, end).';                    % residual
end
