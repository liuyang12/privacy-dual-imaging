%FIG03_TOUCH_OPEN_HAND
% Figure 3 | Imaging privacy threats revealing touch detection in front of 
% the screen. First row: Open hand.
% 
% Citation: Yang Liu, Gregory W. Wornell, William T. Freeman, Fredo Durand.
% Imaging Privacy Threats from an Ambient Light Sensor. preprint
% [submitted]. 2023.
clear; clc;
close all;
% [0] environment configuration
addpath(genpath('../packages')); % packages
addpath('../utils'); % utilities

rawdata_dir = '../rawdata/'; % raw data directory
dataset_dir = '../dataset/'; % dataset directory
results_dir = '../results/'; % results

scfile      = [12]; % scene file(s), average if multiple
bgfile      = [];% background file(s), average if multiple

% [1] parameters
SAMPLE    = 'occluder-mannequin'; % sample (scene) name
h         = 32; % horizontal number of pixels
v         = 32; % vertical number of pixels

distance  = 22; % distance between the sample and the screen [cm]
framerate = 20; % frame rate of the acquisition (camera as single-pixel detector)

BINNING   = 4; % binning size 
deltaH    = 0; % horizontal center frequency (pixel) of the Fourier spectrum
deltaV    = 0; % vertical center frequency (pixel) of the Fourier spectrum
SHAPE     = 'square'; % sampled shape of the Fourier spectrum -- square, circle, or diamond
SHIFTS    = 'four-step'; % four-step or three-step phase shifting for Fourier spectrum acquisition
TYPE      = 'analog'; % type of the pattern -- analog, grayscale, or binary
SAMPRATE  = 1; % sampling rate 
SYMMETRY  = 1; % Hermitian symmetry of real signal
COLORPROJ = false; % color patterns on the screen (projector)
  colormtx = [1 1 0; 0 1 1; 1 0 1]; % RG-GB-BR (Y-C-M) color sampling
  % colormtx = [1 0 0; 0 1 0; 0 0 1]; % R-G-B color sampling
  nch = size(colormtx,1); % number of color channels in color patterns
nColor    = 1; % number of color channels in the measurements
DIFFMEAS  = true; % differential measurement M=(1+M)/2-(1-M)/2
HFLIP     = true; % horizontal flip
UPSIDEDOWN = false; % upside down
RERAWREAD  = false; % re-read the raw measurement data
BRCORRFLAG = true; % brightness correction flag
FIGFLAG    = false; % flag of plotting figure(s)
ALSDATA    = true; % using the real ambient light sensor (ALS) data

% params.BLANKSEG = false; % measurement segmentation using blank frames (sharp drops)
% %   params.nmeas = (v*(h/2+1))*4; % number of measurements (Fourier bases)
%   params.nmeas = v*h*2; % number of measurements (Walsh-Hadamard bases)
  % params.nmeas = 120*2; % 16x16x0.5
%   params.nmeas = 496*2; % 32x32x0.5
% params.interfactor = 1.5; % interval factor (to determine the segmentation)
params.measthresh = 3; % single channel [ambient light sensor output]
% params.avginterval = [0.2, 0.8];
params.processmeth = 'mean'; % processing method to obtain each measurement [maximum of a period]
% params.processmeth = 'max'; % prtocessing method to obtain each measurement [maximum of a period]
% params.avginterval = [0, 1];
% params.avginterval = [-0.2, 1.2];
%   params.patseg = patseg;
if ~exist([dataset_dir SAMPLE], 'dir')
    mkdir([dataset_dir SAMPLE]);
end
if ~exist([results_dir SAMPLE],'dir')
    mkdir([results_dir SAMPLE]);
end

% avrage over all scene measurements
for i = 1:length(scfile)
    ifile = scfile(i);
    
    vname = sprintf('%s%dx%d_%dcm_%dfps_%03d',SAMPLE,h,v,distance,framerate,ifile);
    datafile = [dataset_dir SAMPLE '/' vname '.mat'];

    if exist(datafile, 'file') && (exist('RERAWREAD','var') && ~RERAWREAD)
        load(datafile,'meas');
    else
        % % [2] read the recorded video to the raw correlated intensity measurements
        % I = readvideo2intensity([rawdata_dir '/video/' vname '.avi']);
        % save([rawdata_dir '/mat/' vname],'I');

        % load([rawdata_dir '/mat/' vname],'I');
        
        datfname = [rawdata_dir '/data/' vname '.dat'];
        csvfname = [rawdata_dir '/csv/' vname '.csv'];
        if exist(datfname, 'file')
            ALSDATA = false; % using the camera data
            Iall = load([rawdata_dir '/data/' vname '.dat']);
        elseif exist(csvfname, 'file')
            ALSDATA = true; % using the real ALS data
            % table = readtable(csvfname,'Format','%f%f%d%s');
            table = readtable(csvfname,'Format','%f%d%s');
            Iall = table2array(table(:,1));
%             % RGBW raw data
%             table = readtable(csvfname,'Format','%f%f%f%f%f%s');
%             Iall = table2array(table(:,1:3))';
        else
            error('raw data %s not found.',vname);
        end
        dlen = length(Iall)/nColor;
        I=reshape(Iall,nColor,[])';

        % % [3] turn the raw measurements to the corresponding correlated
        % measurements
        [measraw patseg] = readvidmeasdata(I,params); 

        if DIFFMEAS % differential measurements
            meas = measraw(:,1:2:end) - measraw(:,2:2:end);
        else
            meas = measraw;
        end
        
        if COLORPROJ % color pattern projection
            meas = reshape(meas,[],nch)';
        end

        % save raw data
        save(datafile,'meas');
    end
    
    if i == 1
        scmeas = meas;
    else
        scmeas = (scmeas*(i-1)+meas)/i;
    end
end

% average over all background measurements, if applicable
if ~exist('bgfile','var') || isempty(bgfile)
    meas = scmeas;
else
    for i = 1:length(bgfile)
        ifile = bgfile(i);

        vname = sprintf('%s%dx%d_%dcm_%dfps_%03d',SAMPLE,h,v,distance,framerate,ifile);
        datafile = [dataset_dir SAMPLE '/' vname '.mat'];

        if exist(datafile, 'file') && (exist('RERAWREAD','var') && ~RERAWREAD)
            load(datafile,'meas');
        else
            % % [2] read the recorded video to the raw correlated intensity measurements
            % I = readvideo2intensity([rawdata_dir '/video/' vname '.avi']);
            % save([rawdata_dir '/mat/' vname],'I');

            % load([rawdata_dir '/mat/' vname],'I');

            datfname = [rawdata_dir '/data/' vname '.dat'];
            csvfname = [rawdata_dir '/csv/' vname '.csv'];
            if exist(datfname, 'file')
                ALSDATA = false; % using the camera data
                Iall = load([rawdata_dir '/data/' vname '.dat']);
            elseif exist(csvfname, 'file')
                ALSDATA = true; % using the real ALS data
                % table = readtable(csvfname,'Format','%f%f%d%s');
                table = readtable(csvfname,'Format','%f%d%s');
                Iall = table2array(table(:,1));
                
%                 % RGBW raw data
%                 table = readtable(csvfname,'Format','%f%f%f%f%f%s');
%                 Iall = table2array(table(:,1:3))';
            else
                error('raw data %s not found.',vname);
            end
            dlen = length(Iall)/nColor;
            I=reshape(Iall,nColor,[])';

            % % [3] turn the raw measurements to the corresponding correlated
            % measurements
            [measraw patseg] = readvidmeasdata(I,params); 

            if DIFFMEAS % differential measurements
                meas = measraw(:,1:2:end) - measraw(:,2:2:end);
            else
                meas = measraw;
            end
            
            if COLORPROJ % color pattern projection
                meas = reshape(meas,[],nch)';
            end

            % save raw data
            save(datafile,'meas');
        end

        if i == 1
            bgmeas = meas;
        else
            bgmeas = (bgmeas*(i-1)+meas)/i;
        end
    end
    
    meas = bgmeas - scmeas; % pinhole = background - pinspeck
end

%% [4] reconstruction of the scene, color channel by color channel
%       the denoising-based approximate message passing (D-AMP) algorithm
%       is employed here
csmethod = 'damp'; % compressive sensing (CS) reconstruction method
% csmethod = 'tv'; % compressive sensing (CS) reconstruction method
maxiter = 30; % maximum number of iterations
denoiser = 'FFDNet'; 
% denoiser = 'BM3D'; 
% denoiser = 'TV'; 

% [4.1] load the sensing matrix 
% load('mask.mat');
% load('mask256.mat');

% meas = sum(meas,1);
% nColor = 1;

alpha = 1;
% alpha = 17;
% alpha = 4.2e3;
% alpha = 3.5e3;
% alpha = 2.5e3;
% alpha = 920; % scale coefficient [32x32x100%] occluder-mit [subtraction]
% alpha = 1.21e3; % scale coefficient [32x32x100%] occluder-stripe [subtraction]
% alpha = 1.35e3; % scale coefficient [32x32x100%] occluder-mit [subtraction]
% alpha = 9.33e3; % scale coefficient [32x32x100%] occluder-stripe
% alpha = 2.93e3; % scale coefficient [64x64x50%] occluder-stripe
% alpha = 9.787e4; % scale coefficient [32x32x100%]
% alpha = 3.073e4; % scale coefficient [64x64x50%]
% alpha = 7.521e3; % scale coefficient [128x128x10%]
% alpha = 1.89e3; % scale coefficient [256x256x5%]
% y = meas';
y = meas'/alpha;
% y = sum(meas',2);
% y = (y-mean(y(:)))/alpha;

MAXB = 255;
chnum = size(y,2);

para.rows = v; % [vertical] height
para.cols = h; % [horizontal] width

samprate = round(length(meas)/(h*v)/0.25)*0.25; % sampling rate 

load(sprintf('../packages/data/sim/static/peppers%d_hadamard_samprate%.2f_snr50db.mat',h,samprate),'csdata');
% load(sprintf('../packages/data/sim/static/peppers%d_hadamard_samprate0.50_snr50db.mat',h),'csdata');
% load(sprintf('../packages/data/sim/static/mit_logo%d_video_samprate%.2f_snr50db.mat',h,samprate),'csdata');

% load(sprintf('%s/%s%d_samprate%.2f_snr%ddb.mat',simdatadir,...
%     sparams.sampname,sparams.rows,sparams.samprate,sparams.noisesnr));

%% reconstruction

% csmethod = 'zero-filling';
csmethod = 'pnp';
% csmethod = 'ap';

switch lower(csmethod)
    case 'ap'
        para.DIFFSENS = false; % differential sensing matrix (row subtraction) and the corresponding measurements
        para.csmethod = 'AP'; % alternating projection (AP) method
          para.maxiter = 10; % maximum number of iterations in AP

    case 'zero-filling'
        para.csmethod = 'zero-filling'; % naive zero-filling (not CS though)
                                        %  works only with orthogonal bases
        para.sensmethod = 'hadamard';
        para.sensind = csdata.sensind;
        para.sensmtx = csdata.sensmtx;

    case 'pnp'
        para.csmethod = 'zero-filling'; % naive zero-filling (not CS though)
                                        %  works only with orthogonal bases
        para.sensmethod = 'hadamard';
        para.sensind = csdata.sensind;
        para.sensmtx = csdata.sensmtx;
        % x_zf = cs(csdata.sensmat, y, para);
        [x_zf,Tx_zf] = zero_filling(y, para);
        
        
        para.x0 = x_zf; % zero-filling as initialization
        
        % para.csmethod = 'DAMP'; % D-AMP method (PnP-AMP)
        %   para.denoiser = 'FFDNet'; % FFDNet image denoising
        %   % para.denoiser = 'BM3D'; % FFDNet image denoising
        %   para.maxiter = 5; % maximum number of iterations
        
        para.orthogonal_basis = true;
        
        para.channel_wise = true; % channel-wise reconstruction
        % para.csmethod = 'PnP-ADMM'; % PnP-ADMM method
        %   para.rho = 3.5e2; % multiplier (noise regularier) [32x32x100%] (als) [0,18] lux
        para.csmethod = 'PnP-QCS'; % PnP-ADMM method
          para.p = 3; % norm moment of BPDQ_p
          para.rho = 12e2; % multiplier (noise regularier) [32x32x100%] (als) [0,18] lux
        
        % para.csmethod = 'PnP-GAP'; % PnP-GAP method
        % para.acc = false; % enable accelerated GAP
        
        %   para.denoiser = 'TV'; % total variation denoising
        %   para.maxiter =  50; % maximum number of iterations
          % para.maxiter = 100; % maximum number of iterations
        
        %   para.denoiser = 'BM3D'; % BM3D image denoising
        %   % para.denoiser = 'CBM3D'; % color BM3D image denoising
        %   % para.sigma   = [12  6]/MAXB; % noise deviation (to be estimated and adapted)
        %   % para.maxiter = [10 10];
        %   para.sigma   = [ 6]/MAXB; % noise deviation (to be estimated and adapted)
        %   para.maxiter = [10];
          
        %   para.denoiser = 'WNNM'; % WNNM image denoising
        %   para.range   = 1; % signal value range during recovery
        %   para.sigma   = [12  6]/MAXB; % noise deviation (to be estimated and adapted)
        %   para.maxiter = [20 20];
        
        para.denoiser = 'FFDNet'; % FFDNet image denoising
        if chnum == 1 || (isfield(para,'channel_wise') && para.channel_wise)  % grayscale image
          load(fullfile('models','FFDNet_gray.mat'),'net');
        else % color image
          load(fullfile('models','FFDNet_color.mat'),'net');
        end
        para.net = vl_simplenn_tidy(net);
        para.useGPU = true;
        if para.useGPU
          para.net = vl_simplenn_move(para.net, 'gpu') ;
        end
        para.ffdnetnorm_init = true; % use normalized video for the first 10 iterations
        para.ffdnetnorm = false; % normalize the video before FFDNet video denoising
        para.ffdnetnorm = true; % normalize the video before FFDNet video denoising
        % para.sigma   = [50 25 12  6]/MAXB; % noise deviation (to be estimated and adapted)
        % para.maxiter = [10 10 10 10];
        % para.sigma   = [12  6]/MAXB; % noise deviation (to be estimated and adapted)
        % para.maxiter = [10 10];
        para.sigma   = [12  6]/MAXB; % noise deviation (to be estimated and adapted)
        para.maxiter = [10 10];
    otherwise
        disp('unsupported csmethod.');
end

tstart = tic;
% x = cs(M, y, para); % apply cs method
sensmat = csdata.sensmat;
nsamp = size(y,1);
% nsamp = size(y,1)/8;
x = cs(sensmat(1:nsamp,:), y(1:nsamp,:), para); % apply cs method [integrated channels]
t_cs = toc(tstart);

% if COLORPROJ % color projection
%     x = x*inv(colormtx)';
% end

imraw = reshape(x, [v h chnum]);

% apply intensity compensation to the recovered image
if BRCORRFLAG
    a = 2;
    % b = 0.7;
    b = 1;
    d = h/2.5;
    x0 = round(h/2.2);
    y0 = round(v/3.5);

    x_ = 1:h;
    y_ = 1:v;

    [X,Y] = meshgrid(x_,y_);

    % % linear
    % imraw = imraw.*(1+(abs(X-x0)/a+abs(Y-y0)/b)/d);
    % quadritic
    imraw = imraw.*(1+((X-x0).^2/a^2+(Y-y0).^2/b^2)/d^2);
end

% nColor = min(nColor,size(y,2));
% for icolor = 1:nColor % iterate three color channels
%     yc = y(:,icolor);
%     if ~DIFFMEAS
%         yc = (yc-mean(yc))/a;
%     end
%     switch lower(csmethod)
%         case 'damp'
%             % imraw(:,:,icolor) = DAMP(yc,maxiter,v,h,denoiser,M,[],[]);
%             x = DAMP(yc,maxiter,v,h,denoiser,M,[],[]);
%             imraw(:,:,icolor) = reshape(x, [v h]);
%         case 'tv'
%             opt.row = v;
%             opt.col = h;
%             opt.mu_0 = 1;
%             opt.tol = 1e-3;
%             opt.min_iter = 50;
%             imraw(:,:,icolor) = tvcs(M, yc, opt);
%         otherwise
%             disp('unsupported denoiser!');
%     end
% end

% [5] show the reconstructed scene
% magsize = 16; % magnification factor
magsize = [256 256]; % magnified size

imrecon = imnorm(imraw);
if HFLIP % flip horizontally
    imrecon = flip(imrecon,2);
end
if UPSIDEDOWN % flip horizontally
    imrecon = rot90(imrecon,2);
end
if nColor > 1 || (COLORPROJ && nch > 1)
    blank = zeros(size(imrecon,[1 2]));
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

%% save the reconstruction

if exist('bgfile', 'var') && ~isempty(bgfile)
    if length(scfile) > 1 && length(bgfile) > 1
        vname = sprintf('%s%dx%d_%dcm_%dfps_%03d-%03d_subtracted_by_%03d-%03d',SAMPLE,h,v,distance,framerate,scfile(1),scfile(end),bgfile(1),bgfile(end));
    elseif length(scfile) > 1
        vname = sprintf('%s%dx%d_%dcm_%dfps_ %03d-%03d_subtracted_by_%03d',SAMPLE,h,v,distance,framerate,scfile(1),scfile(end),bgfile(1));
    elseif length(bgfile) > 1
        vname = sprintf('%s%dx%d_%dcm_%dfps_%03d_subtracted_by_%03d-%03d',SAMPLE,h,v,distance,framerate,scfile(1),bgfile(1),bgfile(end));
    else
        vname = sprintf('%s%dx%d_%dcm_%dfps_%03d_subtracted_by_%03d',SAMPLE,h,v,distance,framerate,scfile(1),bgfile(1));
    end
else
    if length(scfile) > 1
        vname = sprintf('%s%dx%d_%dcm_%dfps_%03d-%03d',SAMPLE,h,v,distance,framerate,scfile(1),scfile(end));
    else
        vname = sprintf('%s%dx%d_%dcm_%dfps_%03d',SAMPLE,h,v,distance,framerate,scfile(1));
    end
end

if ALSDATA % using the real ALS data
    vname = [vname '_als'];
else % using camera data 
    vname = [vname '_cam'];
end

if nColor > 1 || (COLORPROJ && nch > 1)
    if isfield(para, 'denoiser')
        if BRCORRFLAG
            savename = sprintf('%s%s/%s_%s-%s_color_brcorr', results_dir, SAMPLE, vname, para.csmethod, para.denoiser);
        else
            savename = sprintf('%s%s/%s_%s-%s_color', results_dir, SAMPLE, vname, para.csmethod, para.denoiser);
        end
    else
        if BRCORRFLAG
            savename = sprintf('%s%s/%s_%s_color_brcorr', results_dir, SAMPLE, vname, para.csmethod);
        else
            savename = sprintf('%s%s/%s_%s_color', results_dir, SAMPLE, vname, para.csmethod);
        end
    end
else
    if isfield(para, 'denoiser')
        if BRCORRFLAG
            savename = sprintf('%s%s/%s_%s-%s_brcorr', results_dir, SAMPLE, vname, para.csmethod, para.denoiser);
        else
            savename = sprintf('%s%s/%s_%s-%s', results_dir, SAMPLE, vname, para.csmethod, para.denoiser);
        end
    else
        if BRCORRFLAG
            savename = sprintf('%s%s/%s_%s_brcorr', results_dir, SAMPLE, vname, para.csmethod);
        else
            savename = sprintf('%s%s/%s_%s', results_dir, SAMPLE, vname, para.csmethod);
        end
    end
end
imwrite(imresize(imrecon,magsize,'nearest'),[savename '.png']);
    saveas(fig,[savename '.fig']);

    
