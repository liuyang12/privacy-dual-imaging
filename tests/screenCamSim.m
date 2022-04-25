%SCREENCAMSIM Simulation of the Screen Camera (ScreeCam).
%   See also GENSCREENCAMPATTERN.

%   Copyright(C) 2019 <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>
%   Last modified Oct 19, 2019.

clear; clc;
% close all;
% [0] environment configuration
addpath('../utils'); % utilities

% rawdata_dir = '../rawdata/'; % rawdata directory
dataset_dir = '../dataset/'; % dataset directory
scene_dir = '../dataset/scenes/static'; % scenes directory
results_dir = '../results/'; % results

% [1] parameters
SAMPLE    = 'map_asia'; % sample used for simulation
h         = 32; % horizontal number of pixels
v         = 32; % vertical number of pixels
BINNING   = 4; % binning size 
deltaH    = 0; % horizontal center frequency (pixel) of the Fourier spectrum
deltaV    = 0; % vertical center frequency (pixel) of the Fourier spectrum
SHAPE     = 'square'; % sampled shape of the Fourier spectrum -- square, circle, or diamond
SHIFTS    = 'four-step'; % four-step or three-step phase shifting for Fourier spectrum acquisition
TYPE      = 'grayscale'; % type of the pattern -- analog, grayscale, or binary
SAMPRATE  = 1; % sampling rate 
SYMMETRY  = 1; % Hermitian symmetry of real signal

BITDEPTH  = 8; % bit depth of the sinusoidal patterns
opts = [];
  opts.method      = 'dither'; % Fourier filtering method
  opts.r_filt      = 1; % pixel radius of the Fourier filter
  opts.cent_includ = 1; % center frequncy included or not. 1 - included, 0 - not included.
  opts.bin_dither  = 0; % convert grayscale to binary using dithering or not. 1 - dithering, 0 - not dithering.
  
% Gaussian blur (mimicking the incoherent light transform process)
GAUSSBLUR = false; 
  SIGMA = 4; % Gaussian blur kernel (sigma)

QUANT = false;
  QUANT_STEPS = 1; % quantilization of the steps of measured values

% [2.1] read the image as the scene
im = imreadallfmt([scene_dir '/' SAMPLE]);
if size(im,3) > 1
    im = rgb2gray(im);
end
im = double(im);
im = imnorm(imresize(im, [v h]));
oim = imnorm(imresize(im,BINNING,'nearest'));

maxB = 2^BITDEPTH-1;

% [2.2] get the single-pixel output from each transform basis
[X,Y] = meshgrid(1:h,1:v);
hc = floor(h/2)+1; % horizontal zero-frequency point after fftshift
vc = floor(v/2)+1; % vertical zero-frequency point after fftshift

VMAX = sum(sum(abs(oim))); % maximum value of the measured signal
fim = zeros(v,h); % Fourier plane
I = zeros(1,4);
% F = zeros(v*(h/2+1)*4,1);
m = 1;
switch SHAPE
    case 'square' % square Fourier spectrum acquisition
        if ~SYMMETRY % full acquisition without using Hermitian symmetry
            r = sqrt(SAMPRATE);
            for jk = -floor(floor(h/2)*r)+deltaH:floor(floor((h-1)/2)*r)+deltaH
                for ik = -floor(floor(v/2)*r)+deltaV:floor(floor((v-1)/2)*r)+deltaV
                    for kphi = 0:1:3
                        switch TYPE % type of patterns
                            case 'analog' % analog patterns
                                P = 1/2+1/2*cos(2*pi/h*jk*X+2*pi/v*ik*Y+kphi*pi/2); % analog
                                if BINNING>1
                                    P = imresize(P,BINNING,'nearest');
                                end
                            case {'gray','grayscale','grey','greyscale'} % grayscale patterns with limited bit depth
                                P = round((1/2+1/2*cos(2*pi/h*jk*X+2*pi/v*ik*Y+kphi*pi/2))*maxB)/maxB; % B-bit grayscale
                                if BINNING>1
                                    P = imresize(P,BINNING,'nearest');
                                end
                            case {'bin','binary'} % binary patterns
                                P = bin2grayfringe(jk,ik,kphi*pi/2,h,v,BINNING,opts); % binary to grascale fringe
                            otherwise
                                error('Unsupported type of patterns -- %s.',TYPE);
                        end
                        if GAUSSBLUR % apply Gaussian blur to the pattern (mimicking the incoherent light transport process)
                            P = imgaussfilt(P, SIGMA);
                        end
                        I(kphi+1) = sum(sum(abs(oim.*P))); % intensity detection (incoherent illumination)
                        F(m) = sum(sum(abs(oim.*P))); % intensity detection (incoherent illumination)
                        if QUANT
                            I(kphi+1) = round(I(kphi+1)/VMAX*QUANT_STEPS)/QUANT_STEPS*VMAX; % quantilization
                            F(m) = round(F(m)/VMAX*QUANT_STEPS)/QUANT_STEPS*VMAX; % quantilization
                        end
                        m = m + 1;
                    end
                    fim(vc+ik-deltaV,hc+jk-deltaH) = (I(1)-I(3)) + 1i*(I(2)-I(4));
                end
            end
        else % half acquisition using Hermitian symmetry
            r = sqrt(SAMPRATE);
            for jk = -floor(floor(h/2)*r)+deltaH:deltaH
                for ik = -floor(floor(v/2)*r)+deltaV:floor(floor((v-1)/2)*r)+deltaV
                    for kphi = 0:1:3
                        switch TYPE % type of patterns
                            case 'analog' % analog patterns
                                P = 1/2+1/2*cos(2*pi/h*jk*X+2*pi/v*ik*Y+kphi*pi/2); % analog
                                if BINNING>1
                                    P = imresize(P,BINNING,'nearest');
                                end
                            case {'gray','grayscale','grey','greyscale'} % grayscale patterns with limited bit depth
                                P = round((1/2+1/2*cos(2*pi/h*jk*X+2*pi/v*ik*Y+kphi*pi/2))*maxB)/maxB; % B-bit grayscale
                                if BINNING>1
                                    P = imresize(P,BINNING,'nearest');
                                end
                            case {'bin','binary'} % binary patterns
                                P = bin2grayfringe(jk,ik,kphi*pi/2,h,v,BINNING,opts); % binary to grascale fringe
                            otherwise
                                error('Unsupported type of patterns -- %s.',TYPE);
                        end
                        if GAUSSBLUR % apply Gaussian blur to the pattern (mimicking the incoherent light transport process)
                            P = imgaussfilt(P, SIGMA);
                        end
                        I(kphi+1) = sum(sum(abs(oim.*P))); % intensity detection (incoherent illumination)
                        F(m) = sum(sum(abs(oim.*P))); % intensity detection (incoherent illumination)
                        if QUANT
                            I(kphi+1) = round(I(kphi+1)/VMAX*QUANT_STEPS)/QUANT_STEPS*VMAX; % quantilization
                            F(m) = round(F(m)/VMAX*QUANT_STEPS)/QUANT_STEPS*VMAX; % quantilization
                        end
                        m = m + 1;
                    end
                    fim(vc+ik-deltaV,hc+jk-deltaH) = (I(1)-I(3)) + 1i*(I(2)-I(4));
                    fim(mod(vc-ik-deltaV-1,v)+1,mod(hc-jk-deltaH-1,h)+1) = (I(1)-I(3)) - 1i*(I(2)-I(4));
                end
            end
        end
end

%% [3] reconstruction of the scene
if ~SYMMETRY % full acquisition without using Hermitian symmetry
    FD = F(1:4:end)-F(3:4:end) + 1i*(F(2:4:end)-F(4:4:end));
    fim2 = reshape(FD, [v,h]);
else % half acquisition using Hermitian symmetry
    FD = F(1:4:end)-F(3:4:end) + 1i*(F(2:4:end)-F(4:4:end));
    FD = reshape(FD, [v,floor(h/2)+1]);
    FD2 = conj(FD([1 v:-1:2],floor((h)/2):-1:2));
    FD = [FD FD2];
    fim2 = reshape(FD, [v,h]);
end
    
scim = ifft2(ifftshift(fim2)); % inverse fast Fourier transform
scim = circshift(scim,[-1 -1]);

% [4] show the recovered image and its transform domain
out = imnorm(abs(scim));
ref = imnorm(imresize(oim,[v h],'nearest'));
sprmse = sqrt(immse(uint8(out*255),uint8(ref*255))); % root mean square error
sppsnr = psnr(uint8(out*255),uint8(ref*255));        % peak signal-to-noise ratio
spssim = ssim(uint8(out*255),uint8(ref*255));        % structure similarity
fig = figure; 
subplot(221); imagesc(ref), colormap(gray), colorbar, axis equal, axis off, title('Ground truth');
subplot(222); imagesc(log(abs(fftshift(fft2(ref)))+1)), colormap(gray), colorbar, axis equal, axis off, title('Fourier spectrum');
subplot(223); imagesc(out), colormap(gray), colorbar, axis equal, axis off, title(sprintf('Reconstructed (PSNR=%.2f dB)',sppsnr));
subplot(224); imagesc(log(abs(fim2)+1)), colormap(gray), colorbar, axis equal, axis off, title('Fourier spectrum');

% [5] save the images and figures
if ~exist([results_dir SAMPLE],'dir')
    mkdir([results_dir SAMPLE]);
end

if QUANT
    if GAUSSBLUR
        saveas(fig,sprintf('%s/%s/%s%dx%d_%s_sigma%d_quant%d.fig',results_dir,SAMPLE,SAMPLE,h,v,TYPE,SIGMA,QUANT_STEPS));
    else
        saveas(fig,sprintf('%s/%s/%s%dx%d_%s_sigma0_quant%d.fig',results_dir,SAMPLE,SAMPLE,h,v,TYPE,QUANT_STEPS));
    end
else
    if GAUSSBLUR
        saveas(fig,sprintf('%s/%s/%s%dx%d_%s_sigma%d.fig',results_dir,SAMPLE,SAMPLE,h,v,TYPE,SIGMA));
    else
        saveas(fig,sprintf('%s/%s/%s%dx%d_%s_sigma0.fig',results_dir,SAMPLE,SAMPLE,h,v,TYPE));
    end
end
