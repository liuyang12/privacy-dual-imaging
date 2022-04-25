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

% [0.1] parameters
h         = 32; % horizontal number of pixels
v         = 32; % vertical number of pixels
SCREENSIZE = [1080 1920]; % size of the screen [height width]
% SCREENSIZE = [v h];
DIFFMEAS = true; % differential measurements
FRAMERATE     = 3; % frame rate of the video
STARTBLANKNUM = 24; % number of blank frames at the beginning
ENDBLANKNUM   = 10; % number of blank frames at the end
PATTNUM  = 4; % number of (same) pattern frames per period
BLANKNUM = 2; % number of blank frames per period
sensmethod = 'hadamard'; % sensing method
transformshape = 'square'; % sampled shape in transform domain
% transformshape = 'zigzag'; % sampled shape in transform domain
samprate = 1; % sampling ratio in transform domain

colormtx = [1 1 1]; % grayscale sampling
% colormtx = [1 0 0; 0 1 0; 0 0 1]; % R-G-B color sampling
% colormtx = [1 1 0; 0 1 1; 1 0 1]; % RG-GB-BR (Y-C-M) color sampling

ncolor = size(colormtx,1); % number of (synthetic) color channels

if exist('sensmethod','var')
    if ncolor>1
        vname = sprintf('%s/%s%dx%d_color_csr%.2f_%dp_%dfps_pat%d',pattern_dir,sensmethod,h,v,samprate,SCREENSIZE(1),FRAMERATE,PATTNUM);
    else
        vname = sprintf('%s/%s%dx%d_csr%.2f_%dp_%dfps_pat%d',pattern_dir,sensmethod,h,v,samprate,SCREENSIZE(1),FRAMERATE,PATTNUM);
    end
else
    if DIFFMEAS
        if ncolor>1
            vname = sprintf('%s/pattdiff%dx%d_color_%dp_%dfps_pat%d',pattern_dir,h,v,SCREENSIZE(1),FRAMERATE,PATTNUM);
        else
            vname = sprintf('%s/pattdiff%dx%d_%dp_%dfps_pat%d',pattern_dir,h,v,SCREENSIZE(1),FRAMERATE,PATTNUM);
        end
    else
        if ncolor>1
            vname = sprintf('%s/pattbinary%dx%d_color_%dp_%dfps_pat%d',pattern_dir,h,v,SCREENSIZE(1),FRAMERATE,PATTNUM);
        else
            vname = sprintf('%s/pattbinary%dx%d_%dp_%dfps_pat%d',pattern_dir,h,v,SCREENSIZE(1),FRAMERATE,PATTNUM);
        end
    end
end

vid = VideoWriter(vname);
vid.FrameRate = FRAMERATE;
open(vid);

% % load the saved sensmat .mat file (to keep exactly the same sensing
% % matrix as the measurement process).
% 
% load('mask.mat','M','B'); % 128x128 sampling rate 0.2
% load('mask256.mat','M','B'); % 256x256 sampling rate 0.05

% generate the sensing matrix according the basis type and sampling
% strategy (to re-genrate the sensing matrix with known configurations)
% 
% [0.3] pre-calculation
switch lower(transformshape)
    case 'square' % center square in tramsform domain
        r = sqrt(samprate); 
    case 'zigzag' % center diamond in transform domain
        if samprate <= 1/2
            r = sqrt(2*samprate); 
        else
            r = 2 - sqrt(2*(1-samprate));
        end
    case 'circular' % center circle in transform domain
        if samprate <= pi/4
            r = 2*sqrt(1/pi*samprate); 
        else
            v = 0:1e-3:pi/2;
            x = (v+cos(v))./(1+sin(v));
            theta = interp1(x, v, samprate, 'linear', 'extrap'); % approximate
            r = 1/cos(pi/4-theta/2) + eps;
        end
    case 'mixed' % mixed orthogonal transform based sensing matrx
        %%% DO NOTHING HERE. %%%
    otherwise
        error('Unsupported shape %s in transform domain!',lower(transformshape));
end

% [1] load or generate data for simulation
rows = h;
cols = v;
nyqnum  = rows*cols; % number of Nyquist samping (N)
sampnum = round(samprate*nyqnum); % sampling number (M)

sparams.orthogonal_basis = false; % [indicator] orthogonal basis or not
% [1.2] generate sensing matrix and get measurement vector as dataset
switch lower(sensmethod)
    case 'binary' % binary pseudo-random sensing matrix (A) 0/1
        sensmat = rand(sampnum,nyqnum)>0.5; % binary sensing matrix (A)
    case 'rademacher' % Rademacher distribution random -1/+1 [near-orthogonal]
        sensmat = 2*(rand(sampnum,nyqnum)>0.5)-1;
    case 'gaussian' % Gaussian sensing matrix (A) -/+
        sensmat = randn(sampnum,nyqnum); 
        % sensmat = 1/sampnum*randn(sampnum,nyqnum); % Gaussian sensing matrix (A)
    case {'fourier', ... % Discrete Fourier transform (DFT) sensing matrix [-1/i,1/i] [orthogonal] 'fourier'='dft'
          'dft', ...     % Discrete Fourier transform (DFT) sensing matrix [-1/i,1/i] [orthogonal] 'fourier'='dft'
          'dct', ...     % Discrete cosine transform (DCT) sensing matrix [-1,1] [orthogonal] 'dct'='cosine'
          'cosine', ...  % Discrete cosine transform (DCT) sensing matrix [-1,1] [orthogonal] 'dct'='cosine'
          'haar', ...    % Haar wavelet transform sensing matrix [-1,1] [orthogonal] 
          'walsh', ...   % Walsh-Hadamard sensing method -1/+1 [orthogonal] 'hadamard'='walsh'
          'hadamard'}    % Walsh-Hadamard sensing method -1/+1 [orthogonal] 'hadamard'='walsh'
        sparams.orthogonal_basis = true;
        x = 1:cols;
        y = 1:rows;
        [X,Y] = meshgrid(x,y);
        rad = eye(nyqnum);
        if sum(strcmpi(sensmethod,{'fourier','dft'})) % Fourier zero-frequency at the center
            r = r/2; % full center
            x0 = ceil((cols+1)/2); % center of columns (center after fftshift) 
            y0 = ceil((rows+1)/2); % center of rows (center after fftshift)
            X = abs(X-x0);
            Y = abs(Y-y0);
        end
        switch lower(transformshape)
            case 'square' % center square in tramsform domain
                mask = (X/cols<=r) & (Y/rows<=r); 
            case 'zigzag' % center diamond in transform domain
                mask = X/cols+Y/rows<=r;
            case 'circular' % center circle in transform domain
                mask = (X/cols).^2+(Y/rows).^2<=r^2; 
            case 'mixed' % mixed Walsh-Hadamard sensing matrx
                mask = zeros(rows,cols);
                allrowset = randperm(nyqnum);
                mask(allrowset(randperm(sampnum))) = 1;
                rad = diag((rand(nyqnum,1)>0.5)*2-1);
            otherwise
                error('Unsupported shape %s in transform domain!',lower(sparams.transformshape));
        end
        sensind = find(mask>0); % index of selected rows
        sampnum = length(sensind); % revise number of samples
        switch lower(sensmethod)
            case {'fourier','dft'} % Discrete Fourier transform (DFT) sensing matrix [-1/i,1/i] [orthogonal]
                Wr = dftmtx(rows); % DFT matrix in rows
                Wc = dftmtx(cols); % DFT matrix in columns
                sensmtx = kron(Wc,Wr); % full DFT matrix 
                shiftind = fftshift(reshape(1:rows*cols,[rows,cols])); % shift index of fftshift
                sensmtx = sensmtx(shiftind,:); % fftshift for the rows (zero frequency in the center)
            case {'dct','cosine'} % Discrete cosine transform (DCT) sensing matrix [-1,1] [orthogonal]
                Wr = dctmtx(rows); % DCT matrix in rows
                Wc = dctmtx(cols); % DCT matrix in columns
                % sensmtx = kron(Wc,Wr); % full DCT matrix 
                sensmtx = kron(Wc,Wr)*sqrt(rows*cols)/2; % full DCT matrix [-1,1]
            case {'haar'} % Haar wavelet transform sensing matrix [-1,1] [orthogonal] 
                Wr = haarmtx(rows); % Haar wavelet transform matrix in rows
                Wc = haarmtx(cols); % Haar wavelet transform matrix in columns
                sensmtx = kron(Wc,Wr)*2; % full Haar wavelet transform matrix 
            case {'walsh','hadamard'} % Walsh-Hadamard sensing method -1/+1 [orthogonal]
                Wr = walsh(rows); % Walsh-Hadamard matrix in rows
                Wc = walsh(cols); % Walsh-Hadamard matrix in columns
                sensmtx = kron(Wc,Wr); % full Walsh-Hadamard matrix 
        end
        sensmat = sensmtx(sensind,:)*rad;
        csdata.sensmtx = sensmtx; % save the full sensing matrix for zero-filling
        csdata.sensind = sensind; % save indices of selected bases for zero-filling
    case 'mixhadamard' % mixed Hadamard sensing method -1/+1 [orthogonal]
        sparams.orthogonal_basis = true;
        Hr = hadamard(rows); % Hadamard matrix in rows
        Hc = hadamard(cols); % Hadamard matrix in columns
        Had = kron(Hc',Hr');
        Rad = double(rand(nyqnum,1)>0.5)*2-1;
        allrowset = randperm(nyqnum);
        rowset = allrowset(randperm(sampnum));
        % rowset = randperm(sampnum);
        sensmat = Had(rowset,:)*diag(Rad);
        % sensmat = Had(rowset,:);
        csdata.sensmtx = Had;
        csdata.sensind = rowset;
        csdata.Rad    = Rad;
    case 'video' % generate sensing matrix from a pre-defined video sequence
        vid = VideoReader(sparams.vpath);
        nframe = vid.NumFrames; % number of frames in the video
        w = vid.Width;
        h = vid.Height;
        if ~isfield(sparams, 'startframe'), sparams.startframe = 1; end
        if ~isfield(sparams, 'step'),             sparams.step = 1; end
        if ~isfield(sparams, 'cropsize'),     sparams.cropsize = floor(min(h/rows,w/cols)*[rows cols]); end
        sparams.endframe = sparams.startframe + (sampnum-1)*sparams.step;
        ch = sparams.cropsize(1); cw = sparams.cropsize(2);
        j = 1;
        sensmat = zeros(sampnum, nyqnum);
        % for iframe = sparams.startframe:sparams.step:sparams.endframe
        iframe = sparams.startframe;
        while (j <= sampnum) % && hasFrame(vid)
            img_color = read(vid, iframe);
            img_gray = im2double(rgb2gray(img_color));
            img = img_gray(floor((h-ch)/2)+(1:ch), floor((w-cw)/2)+(1:cw));
            if isfield(sparams,'SIMREALSIZE') && sparams.SIMREALSIZE
                meas_nonoise(j,:) = img(:)'*samp_rzvec;
            end
            % img = imresize(img, [rows, cols], 'bilinear'); % [OPT#1] image resizing
            img = imbinning(img, [ch/rows cw/cols]); % [OPT#2] pixel binning
            % img = round(img*255)/255;
            if mean(img(:)/(ch*cw)*rows*cols) > 0.1 % omit (literally) blank frames
                sensmat(j, :) = img(:);    
                j = j+1;
            else
                fprintf('skipped frame #%d, mean %.4f.\n',iframe, mean(img(:)/(ch*cw)*rows*cols));
            end
            iframe = iframe + sparams.step;
        end
        iframe
        
    otherwise
        error('unsupported sensing method %s.\n', sensmethod);
end

N = size(sensmat,1); % number of patterns in the pile

% [2] write each pattern and blank ones into the video
blank = zeros(SCREENSIZE,'uint8'); % blank frame
newframe = zeros([SCREENSIZE 3],'uint8'); % new frame (consider RGB color)
for ii = 1:STARTBLANKNUM
    writeVideo(vid,blank);
end

for icolor = 1:ncolor
    color = colormtx(icolor,:); % RGB color for each pattern
    for n = 1:N
        B = 1/2*(sensmtx(n,:)+1);
        P = reshape(B,[v h]); % B=(1+M)/2
        sbin = floor(SCREENSIZE(1)/v);
        patt = imresize(uint8(imnorm(P)*255),sbin,'nearest');
        newpatt = zeros([size(patt) 3],'uint8'); % three color channels
        for c = 1:3
            newpatt(:,:,c) = patt*color(c);
        end
        newframe(round((SCREENSIZE(1)-sbin*v)/2)+1:round((SCREENSIZE(1)-sbin*v)/2)+sbin*v,round((SCREENSIZE(2)-sbin*h)/2)+1:round((SCREENSIZE(2)-sbin*h)/2)+sbin*h,:) = newpatt;
        for ii = 1:PATTNUM % write pattern frames
            writeVideo(vid,newframe);
        end
        for ii = 1:BLANKNUM % write blank frames
            writeVideo(vid,blank);
        end
        if DIFFMEAS % differential measurements
            P = reshape(1-B,[h v]); % 1-B=(1-M)/2
            sbin = floor(SCREENSIZE(1)/v);
            patt = imresize(uint8(imnorm(P)*255),sbin,'nearest');
            newpatt = zeros([size(patt) 3],'uint8'); % three color channels
            for c = 1:3
                newpatt(:,:,c) = patt*color(c);
            end
            newframe(round((SCREENSIZE(1)-sbin*v)/2)+1:round((SCREENSIZE(1)-sbin*v)/2)+sbin*v,round((SCREENSIZE(2)-sbin*h)/2)+1:round((SCREENSIZE(2)-sbin*h)/2)+sbin*h,:) = newpatt;
            for ii = 1:PATTNUM % write pattern frames
                writeVideo(vid,newframe);
            end
            for ii = 1:BLANKNUM % write blank frames
                writeVideo(vid,blank);
            end
        end
    end
end

for ii = 1:ENDBLANKNUM
    writeVideo(vid,blank);
end

close(vid); % close video

