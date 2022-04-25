%GENSCREENCAMFILM Generate the patterns for Screen Camera (ScreenCam).
%   See also SCREENCAMSIM.

%   Copyright(C) 2019 <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>
%   Last modified July 8, 2020.

clear; clc;
% close all;
% [0] environment configuration
addpath('../utils'); % utilities

pattern_dir = '../pattern/'; % pattern directory
simdata_dir = '../packages/data/'; % simulation data directory

% [0.1] load simulation preferences
load(sprintf('%s/sim_prefs.mat',simdata_dir));

% [1] parameters
% filmname = 'ready_player_one'; % name of the film
filmname = 'tom_and_jerry'; % name of the film
sensmethod = sparams.sensmethod;  % sensing method
samprate   = sparams.samprate;    % sampling ratio in transform domain
h          = sparams.cropsize(2); % horizontal number of pixels
v          = sparams.cropsize(1); % vertical number of pixels

SCREENSIZE = [1080 1920]; % size of the screen [height width]
% SCREENSIZE = [v h];
DIFFMEAS = false; % differential measurements
% FRAMERATE = 30; % frame rate of the video
FRAMERATE = 12; % frame rate of the video
STARTBLANKNUM = 30; % number of blank frames at the beginning
ENDBLANKNUM   = 15; % number of blank frames at the end
PATTNUM  = 4; % number of (same) pattern frames per period
BLANKNUM = 2; % number of blank frames per period

vname = sprintf('%s/%s%dx%d_csr%.2f_%dp_%dfps_pat%d',pattern_dir,filmname,h,v,samprate,SCREENSIZE(1),FRAMERATE,PATTNUM);

vid = VideoWriter(vname);
vid.FrameRate = FRAMERATE;
open(vid);

rows = sparams.rows; % number of rows in the recovery pixel resolution
cols = sparams.cols; % number of rows in the recovery pixel resolution
N = round(samprate*rows*cols); % number of patterns in the pile

% [1.1] read the film and pre-set parameters
film = VideoReader(sparams.vpath);
W = film.Width;  % width of the film sequence
H = film.Height; % height of the flim sequence

if ~isfield(sparams, 'startframe'), sparams.startframe = 1; end
if ~isfield(sparams, 'step'),             sparams.step = 1; end
if ~isfield(sparams, 'cropsize'),     sparams.cropsize = floor(min(h/rows,w/cols)*[rows cols]); end 
ch = sparams.cropsize(1); cw = sparams.cropsize(2);

% [2] write each pattern and blank ones into the video
blank = zeros(SCREENSIZE,'uint8'); % blank frame
newframe = zeros(SCREENSIZE,'uint8'); % new frame
for ii = 1:STARTBLANKNUM
    writeVideo(vid,blank);
end

n = 1;
iframe = sparams.startframe;
while (n <= N)
    img_color = read(film, iframe);
    img_gray = rgb2gray(img_color);
    img = img_gray(floor((H-ch)/2)+(1:ch), floor((W-cw)/2)+(1:cw));
    
    if mean(img(:)) > 0.1*255
        n = n+1;
        sbin = floor(SCREENSIZE(1)/v);
        patt = imresize(img,sbin,'nearest');
        newframe(round((SCREENSIZE(1)-sbin*v)/2)+1:round((SCREENSIZE(1)-sbin*v)/2)+sbin*v,round((SCREENSIZE(2)-sbin*h)/2)+1:round((SCREENSIZE(2)-sbin*h)/2)+sbin*h) = patt;
        for ii = 1:PATTNUM % write pattern frames
            writeVideo(vid,newframe);
        end
        for ii = 1:BLANKNUM % write blank frames
            writeVideo(vid,blank);
        end
    else
        fprintf('skipped frame #%d, mean %.4f.\n',iframe, mean(img(:)));
    end
    iframe = iframe + sparams.step;
end
iframe

for ii = 1:ENDBLANKNUM
    writeVideo(vid,blank);
end

close(vid); % close video

