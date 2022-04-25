%SCREENCAMRECON Reconstruction of the Screen Camera. 
clear; clc;
% close all;
% [0] environment configuration
addpath('../utils'); % utilities

% rawdata_dir = '../../rawdata/'; % raw data directory
rawdata_dir = '/data/yliu/docs/Dropbox (MIT)/proj/screencam/rawdata/'; % raw data directory
dataset_dir = '../dataset/'; % dataset directory
results_dir = '../results/'; % results

ifile     = 21;

% [1] parameters
% SAMPLE    = 'mit'; % sample (scene) name
% SAMPLE    = 'apple_logo'; % sample (scene) name
% SAMPLE    = 'vase'; % sample (scene) name
% SAMPLE    = 'chipbag'; % sample (scene) name
% SAMPLE    = 'csail'; % sample (scene) name
% SAMPLE    = 'cup'; % sample (scene) name
% SAMPLE    = 'occluder-lego'; % sample (scene) name
% SAMPLE    = 'box'; % sample (scene) name
% SAMPLE    = 'mirror-rect'; % sample (scene) name
% SAMPLE    = 'eagle'; % sample (scene) name
SAMPLE    = 'occluder-stripe'; % sample (scene) name
h         = 32; % horizontal number of pixels
v         = 32; % vertical number of pixels

distance  = 22; % distance between the sample and the screen [cm]
% framerate = 174; % frame rate of the acquisition (camera as single-pixel detector)
framerate = 360; % frame rate of the acquisition (camera as single-pixel detector)


BINNING   = 4; % binning size 
deltaH    = 0; % horizontal center frequency (pixel) of the Fourier spectrum
deltaV    = 0; % vertical center frequency (pixel) of the Fourier spectrum
SHAPE     = 'square'; % sampled shape of the Fourier spectrum -- square, circle, or diamond
SHIFTS    = 'four-step'; % four-step or three-step phase shifting for Fourier spectrum acquisition
TYPE      = 'analog'; % type of the pattern -- analog, grayscale, or binary
SAMPRATE  = 1; % sampling rate 
SYMMETRY  = 1; % Hermitian symmetry of real signal
nColor    = 3; % number of color channels
HFLIP     = true; % horizontal flip

% % [2] read the recorded video to the raw correlated intensity measurements
vname = sprintf('%s%dx%d_%dcm_%dfps_%03d',SAMPLE,h,v,distance,framerate,ifile);

% I = readvideo2intensity([rawdata_dir '/video/' vname '.avi']);
% save([rawdata_dir '/mat/' vname],'I');

% load([rawdata_dir '/mat/' vname],'I');

Iall=load([rawdata_dir '/data/' vname '.dat']);
dlen = length(Iall)/nColor;
I=reshape(Iall,nColor,[])';

%% [3] turn the raw measurements to the corresponding correlated
if ~exist([dataset_dir SAMPLE], 'dir')
    mkdir([dataset_dir SAMPLE]);
end
if ~exist([results_dir SAMPLE],'dir')
    mkdir([results_dir SAMPLE]);
end

COLOR = true;
FIGFLAG = false;
% nColor = 1; % number of color channels

% measurements
params.measthresh = 1e7; % three channels
params.nmeas = (v*(h/2+1))*4; % number of measurements
meas = readvidmeasdata(I,params); 

% load data
save([dataset_dir SAMPLE '/' vname],'meas');

% load([dataset_dir SAMPLE '/' vname],'meas');

%% [4] reconstruction of the scene, color channel by color channel
magsize = 8;

for iColor = 1:nColor
    meas_ch = meas(iColor,:);
    
    
    if ~SYMMETRY % full acquisition without using Hermitian symmetry
        FD = meas_ch(1:4:end)-meas_ch(3:4:end) + 1i*(meas_ch(2:4:end)-meas_ch(4:4:end));
        fim2 = reshape(FD, [v,h]);
    else % half acquisition using Hermitian symmetry
        FD = meas_ch(1:4:end)-meas_ch(3:4:end) + 1i*(meas_ch(2:4:end)-meas_ch(4:4:end));
        FD = reshape(FD, [v,floor(h/2)+1]);
        FD2 = conj(FD([1 v:-1:2],floor((h)/2):-1:2));
        FD = [FD FD2];
        fim2 = reshape(FD, [v,h]);
    end

    scim = ifft2(ifftshift(fim2)); % inverse fast Fourier transform
    scim = circshift(scim,[-1 -1]);

    % [5] show the recovered image and its transform domain
    out = imnorm(abs(scim));
    imrecon(:,:,iColor) = abs(scim);
    % imrecon(:,:,iColor) = out;
    
    % imwrite(imresize(out,magsize,'nearest'),[results_dir SAMPLE '/' vname '.png']);
    
    if FIGFLAG
        fig = figure; 
        subplot(221); imagesc(out), colormap(gray), colorbar, axis equal, axis off, title(sprintf('Reconstructed'));
        subplot(222); imagesc(log(abs(fim2)+1)), colormap(gray), colorbar, axis equal, axis off, title('Fourier spectrum');

        % saveas(fig,[results_dir SAMPLE '/' vname '.fig']);
    end
    
end

%% [5] show the reconstructed scene
magsize = 8;

imrecon = imnorm(imrecon);
if HFLIP % flip horizontally
    imrecon = flip(imrecon,2);
end
% imrecon = imrecon*1.5;
if nColor > 1
    blank = zeros(size(scim));
    red_channel = cat(3, imrecon(:,:,1), blank, blank);
    green_channel = cat(3, blank, imrecon(:,:,2), blank);
    blue_channel = cat(3, blank, blank, imrecon(:,:,3));

    fig = figure;
    subplot(221); imshow(imresize(red_channel,magsize,'nearest')); title('Red channel');
    subplot(222); imshow(imresize(green_channel,magsize,'nearest')); title('Green channel');
    subplot(223); imshow(imresize(blue_channel,magsize,'nearest')); title('Blue channel');
    subplot(224); imshow(imresize(imrecon,magsize,'nearest')); title('RGB color');
else
    fig = figure;
    imshow(imresize(imrecon,magsize,'nearest')); title('Reconstructued');
end

if nColor > 1
    imwrite(imresize(imrecon,magsize,'nearest'),[results_dir SAMPLE '/' vname '_color.png']);
    saveas(fig,[results_dir SAMPLE '/' vname '_color.fig']);
else
    imwrite(imresize(imrecon,magsize,'nearest'),[results_dir SAMPLE '/' vname '.png']);
    saveas(fig,[results_dir SAMPLE '/' vname '.fig']);
end

