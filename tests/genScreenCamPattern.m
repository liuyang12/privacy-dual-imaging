%GENSCREENCAMPATTERN Generate the patterns for Screen Camera (ScreenCam).
%   See also SCREENCAMSIM.

%   Copyright(C) 2019 <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>
%   Last modified Oct 19, 2019.

clear; clc;
% close all;
% [0] environment configuration
addpath('../utils'); % utilities

pattern_dir = '../pattern/'; % pattern directory
dataset_dir = '../dataset/'; % dataset directory
scene_dir = '../dataset/scenes/static'; % scenes directory

% [1] parameters
SAMPLE    = 'lena'; % sample used for simulation
h         = 64; % horizontal number of pixels
v         = 64; % vertical number of pixels
BINNING   = 1; % binning size 
deltaH    = 0; % horizontal center frequency (pixel) of the Fourier spectrum
deltaV    = 0; % vertical center frequency (pixel) of the Fourier spectrum
SHAPE     = 'square'; % sampled shape of the Fourier spectrum -- square, circle, or diamond
SHIFTS    = 'four-step'; % four-step or three-step phase shifting for Fourier spectrum acquisition
TYPE      = 'analog'; % type of the pattern -- analog, grayscale, or binary
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
  SIGMA = 5; % Gaussian blur kernel (sigma)

SCREENSIZE = [1080 1920];
FRAMERATE = 30; % frame rate of the video
STARTBLANKNUM = 60; % number of blank frames at the beginning
ENDBLANKNUM = 10; % number of blank frames at the end
BLANKNUM = 1; % number of blank frames per period
PATTNUM = 2; % number of (same) pattern frames per period

vid = VideoWriter(sprintf('%s/patt%dx%d_%dp_%dfps_pat%d',pattern_dir,h,v,SCREENSIZE(1),FRAMERATE,PATTNUM));
vid.FrameRate = FRAMERATE;
open(vid);

% [2.1] read the image as the scene
im = imreadallfmt([scene_dir '/' SAMPLE]);
im = double(im);
im = imnorm(imresize(im, [v h]));
oim = imnorm(imresize(im,BINNING,'nearest'));

maxB = 2^BITDEPTH-1;

% [2.2] get the single-pixel output from each transform basis
[X,Y] = meshgrid(1:h,1:v);
hc = floor(h/2)+1; % horizontal zero-frequency point after fftshift
vc = floor(v/2)+1; % vertical zero-frequency point after fftshift

fim = zeros(v,h); % Fourier plane
I = zeros(1,4);
% F = zeros(v*(h/2+1)*4,1);
m = 1;
kp = 1;
blank = zeros(SCREENSIZE,'uint8'); % blank frame
newframe = zeros(SCREENSIZE,'uint8'); % new frame
for ii = 1:STARTBLANKNUM
    writeVideo(vid,blank);
end
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
                        sbin = floor(SCREENSIZE(1)/v);
                        patt = imresize(uint8(imnorm(P)*255),sbin,'nearest');
                        newframe(round((SCREENSIZE(1)-sbin*v)/2):round((SCREENSIZE(1)-sbin*v)/2)+sbin*v-1,round((SCREENSIZE(2)-sbin*h)/2):round((SCREENSIZE(2)-sbin*h)/2)+sbin*h-1) = patt;
                        for ii = 1:PATTNUM % write pattern frames
                            writeVideo(vid,newframe);
                        end
                        for ii = 1:BLANKNUM % write blank frames
                            writeVideo(vid,blank);
                        end
                        
                        if GAUSSBLUR % apply Gaussian blur to the pattern (mimicking the incoherent light transport process)
                            P = imgaussfilt(P, SIGMA);
                        end
                        I(kphi+1) = sum(sum(abs(oim.*P))); % intensity detection (incoherent illumination)
                        F(m) = sum(sum(abs(oim.*P))); % intensity detection (incoherent illumination)
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
                        
                        sbin = floor(SCREENSIZE(1)/v);
                        patt = imresize(uint8(imnorm(P)*255),sbin,'nearest');
                        newframe(round((SCREENSIZE(1)-sbin*v)/2):round((SCREENSIZE(1)-sbin*v)/2)+sbin*v-1,round((SCREENSIZE(2)-sbin*h)/2):round((SCREENSIZE(2)-sbin*h)/2)+sbin*h-1) = patt;
                        for ii = 1:PATTNUM % write pattern frames
                            writeVideo(vid,newframe);
                        end
                        for ii = 1:BLANKNUM % write blank frames
                            writeVideo(vid,blank);
                        end
                        
                        if GAUSSBLUR % apply Gaussian blur to the pattern (mimicking the incoherent light transport process)
                            P = imgaussfilt(P, SIGMA);
                        end
                        I(kphi+1) = sum(sum(abs(oim.*P))); % intensity detection (incoherent illumination)
                        F(m) = sum(sum(abs(oim.*P))); % intensity detection (incoherent illumination)
                        m = m + 1;
                    end
                    fim(vc+ik-deltaV,hc+jk-deltaH) = (I(1)-I(3)) + 1i*(I(2)-I(4));
                    fim(mod(vc-ik-deltaV-1,v)+1,mod(hc-jk-deltaH-1,h)+1) = (I(1)-I(3)) - 1i*(I(2)-I(4));
                end
            end
        end
end

for ii = 1:ENDBLANKNUM
    writeVideo(vid,blank);
end

close(vid); % close video

