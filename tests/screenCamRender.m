
clear; clc; 
% close all

imdir = '../../sim/results';

bwlist = {'black','white'};
imno = 26;
bw   =  2; % black or white scene background

magsize = [256 256];
imI = imread(sprintf('%s/screencam_occluder_i_%s_%02d.tif',imdir,bwlist{bw},imno));
imT = imread(sprintf('%s/screencam_occluder_t_%s_%02d.tif',imdir,bwlist{bw},imno));

im = flipud(double(imI - imT));

figure; imshow(im/max(im(:)));

imwrite(im/max(im(:)),sprintf('%s/screencam_subtracted_%s_%02d.tif',imdir,bwlist{bw},imno));
