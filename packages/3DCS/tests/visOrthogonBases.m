%VISORTHOGONBASES Visualize Orthogonal Bases. 
%   See also TESTCS, GENDATA, TIGHT_SUBPLOT.

% clear; clc;
% close all;

% [1] read the sensing matrix
% sensmat

basisRange = [8 8]; % [nrows ncols] to be displayed
rRange = basisRange(1);
cRange = basisRange(2);

% [2] visualize each row of the sensing matrix and arrange them as a square
fig = figure('position', [100 100 800 800]);
hplot = tight_subplot(rRange,cRange,[0.01 0.01], [0.01, 0.01], [0.01 0.01]);

for irow = 1:rRange
    for icol = 1:cRange
        vect = sensmat((irow-1)*cols+icol,:);
        if isreal(sensmat) % real-valued sensing matrix [-1,1]
            if min(vect) >= 0 % [0,1]
                vectnorm = vect;
            else % [-1,1]
                vectnorm = (1+vect)/2; % matlab vec() order row-first
            end
        else % complex-valued sensing matrix [-1/-i,1/i] 
            vectnorm = (1+cos(angle(vect)))/2; % matlab vec() order row-first
        end
        patt = reshape(vectnorm, [rows,cols]);
        axes(hplot(irow+(icol-1)*rRange)); % subplot order column-first
        imshow(patt);
    end
end

% [3] save the figure to desired directory
outdir = '../out';  % simulation output directory
figdir = [outdir '/fig']; % figure directory
if ~exist(figdir, 'dir')
    mkdir(figdir);
end

saveas(fig, sprintf('%s/%s%dx%d_visual.svg',...
                     figdir,sparams.sensmethod,cols,rows));

