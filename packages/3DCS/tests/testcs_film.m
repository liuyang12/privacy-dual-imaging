%TESTCS_FILM test compressive sensing (CS) algorithms for sparse 
%reconstruction, where the sensing matrix is formed by a pre-defined film
%or video sequence.
%   See also CS, TEST_CS.
%% step 1. generate and save simulation dataset
clear; 
% close all; clc;
% [1.0] parameters for simulation
sparams = []; % paramters for simulation
    sparams.rows = 32; % number of rows (m)
    sparams.cols = 32; % number of columns (n)
    sparams.samprate = 8; % sampling rate (gamma)
    sparams.noisesnr = 50; % signal-to-noise ratio (sigma)
    % sparams.noisesnr = inf; % signal-to-noise ratio (sigma)
    % sparams.sampname = 'gi'; % binary sample 
    % sparams.sampname = 'mit_logo'; % binary sample 
    % sparams.sampname = 'mit_logo_sim'; % binary sample 
    % sparams.sampname = 'hand_touch'; % binary sample 
    sparams.sampname = 'touch_real'; % binary sample 
    % sparams.sampname = 'apple_logo'; % binary sample 
    % sparams.sampname = 'eagle'; % binary sample 
    % sparams.sampname = 'angrybird'; % gray-scale sample 
    % sparams.sampname = 'phantom'; % simple gray-scale sample
    % sparams.sampname = 'cameraman'; % simple gray-scale sample
    % sparams.sampname = 'peppers'; % simple gray-scale sample
    % sparams.sensmethod = 'binary'; % sensing matrix (multiple patterns) is binary
    % sparams.sensmethod = 'gaussian'; % sensing matrix (multiple patterns) is -/+
    % sparams.sensmethod = 'rademacher'; % sensing matrix (multiple patterns) is -1/+1
    % sparams.sensmethod = 'hadamard'; % sensing matrix (multiple patterns) is -1/+1
    % sparams.transformshape = 'square'; % sampled shape in transform domain
    % sparams.sensmethod = 'mixhadamard'; % sensing matrix (multiple patterns) is -1/+1
    sparams.sensmethod = 'video'; % sensing matrix from video [0,1]
    % sparams.vpath = '/proj/screencam/video/ready_player_one_last_fight_240p.mp4'; % video sequence as sensing matrix
    % sparams.startframe = 400; % start frame number
    % sparams.cropsize = [160 160]; % crop size of the original video
    sparams.vpath = '/proj/screencam/video/tom_and_jerry_chase_240p.mp4'; % video sequence as sensing matrix
    sparams.filmname = 'tom_and_jerry';
    sparams.startframe = 1; % start frame number
    sparams.cropsize = [160 160]; % crop size of the original video
    sparams.step = 1; % step of the video sequence in forming sensing matrix

    % sparams.vpath = '/proj/screencam/video/fireworks_240p.mp4'; % video sequence as sensing matrix
    % sparams.filmname = 'fireworks';
    % sparams.startframe = 1; % start frame number
    % sparams.cropsize = [160 160]; % crop size of the original video
    % sparams.step = 2; % step of the video sequence in forming sensing matrix

    % sparams.vpath = '/proj/screencam/video/fpv_fireworks_240p.mp4'; % video sequence as sensing matrix
    % sparams.filmname = 'fpv_fireworks';
    % sparams.startframe = 2525; % start frame number
    % sparams.cropsize = [160 160]; % crop size of the original video
    % sparams.step = 1; % step of the video sequence in forming sensing matrix

    % sparams.vpath = '/proj/screencam/video/fpv_racing_240p.mp4'; % video sequence as sensing matrix
    % sparams.filmname = 'fpv_racing';
    % sparams.startframe = 550; % start frame number
    % sparams.cropsize = [160 160]; % crop size of the original video
    % sparams.step = 1; % step of the video sequence in forming sensing matrix

    % sparams.vpath = '/proj/screencam/video/fpv_mountain_240p.mp4'; % video sequence as sensing matrix
    % sparams.filmname = 'fpv_mountain';
    % % sparams.startframe = 4000; % start frame number - start from 4000  x1 [9.23 dB]
    % sparams.startframe = 100; % start frame number - start from 4000 x8 [9.23 dB]
    % sparams.cropsize = [160 160]; % crop size of the original video
    % sparams.step = 2; % step of the video sequence in forming sensing matrix

    % sparams.vpath = '/proj/screencam/video/speedriding_240p.mp4'; % video sequence as sensing matrix
    % sparams.filmname = 'speedriding';
    % sparams.startframe = 100; % start frame number - start from 4000 x8 [9.23 dB]
    % sparams.cropsize = [160 160]; % crop size of the original video
    % sparams.step = 1; % step of the video sequence in forming sensing matrix

    sparams.SIMREALSIZE    = true;  % use the real size as the cropsize

    % sparams.SAMPLE_RGB     = false; % sample is binary (1-Y, 0-N)
    sparams.SAMPLE_BINARY  = false; % sample is binary (1-Y, 0-N)
    % sparams.savesize = 2; % saving size of the image (magnification)
    sparams.savesize = [256,256]; % saving size of the image (size)
    % sparams.binsize = 1; % binning size of the spatial light modulator
    sparams.QUANTIZATION = true; % apply quantization to the direct measurement
    sparams.QUANT_STEPS  = 63; % steps of quantization

simdatadir = '../../data/sim/static';   % simulation data directory
if strcmpi(sparams.sensmethod, 'video') && isfield(sparams, 'filmname')
    simfile = sprintf('%s/%s%d_%s_%s_samprate%.2f_snr%ddb.mat',simdatadir,sparams.sampname,sparams.rows,sparams.sensmethod,sparams.filmname,sparams.samprate,sparams.noisesnr);
else
    simfile = sprintf('%s/%s%d_%s_samprate%.2f_snr%ddb.mat',simdatadir,sparams.sampname,sparams.rows,sparams.sensmethod,sparams.samprate,sparams.noisesnr);
end

REGENDATA = true; % flag of regenerating date (1-Y,0-N)
% REGENDATA = false; % flag of regenerating date (1-Y,0-N)
if REGENDATA || ~exist(simfile, "file")
% [1] generate and save dataset
    gendata(sparams);
end

%% step 2. load and apply CS method to the dataset
% [2] load dataset
% clear; 
% clc;
% close all;
addpath('../utils')
addpath('../algorithms')
addpath(genpath('../packages'))

datadir    = '../../data';       % data directory

outdir     = '../out';           % simulation output directory
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% [2.1] load simulation preferences
load(sprintf('%s/sim_prefs.mat',datadir));
% [2.2] load dataset
load(simfile)

samp    = csdata.samp;    % sample
sensmat = csdata.sensmat; % sensing matrixs
meas    = csdata.meas;    % measurement vector

rows    = sparams.rows;
cols    = sparams.cols;
MAXB    = 255;           % maximum pixel value of the image (8-bit -> 255)
chnum   = size(meas,2);  % number of color channels (1 for grayscale, 3 for color RGB)

fprintf('Condition number of sensing matrix: %.1f\n', cond(double(sensmat)));
%% [3] apply CS method
csparams = []; % parameters for CS method
    csparams.rows     = rows;
    csparams.cols     = cols;
    % % csparams.orthogonal_basis = sparams.orthogonal_basis;
    % csparams.use_eigendecomp = true;
    %   csparams.condnumb_upper = 2e2;
    %   csparams.orthogonal_basis = true;
    csparams.use_eigendecomp = false;
      csparams.orthogonal_basis = false;
    % csparams.DIFFSENS = true; % differential sensing matrix (row subtraction) and the corresponding measurment
    % csparams.srbasis  = 'haar'; % sparse representation basis
    % csparams.srbasis  = 'daubechies'; % sparse representation basis
    % csparams.srbasis  = 'dct'; % sparse representation basis

    % csparams.csmethod = 'GPSR'; % GPSR ell_1 solver
    % csparams.csmethod = 'TGI'; % traditional correlation-based ghost imaging method
    csparams.csmethod = 'AP'; % alternating projection (AP) method
      csparams.maxiter = 100; % maximum number of iterations in AP
    % csparams.csmethod = 'TV'; % TV regularization method
    % csparams.csmethod = 'GAP'; % GAP ell_1 solver
    % csparams.csmethod = 'GAP-TV'; % GAP-TV  solver
    % csparams.csmethod = 'GAP2D'; % GAP2D solver

%     csparams.csmethod = 'zero-filling'; % naive zero-filling (not CS though)
%                                         %  works only with orthogonal bases
%       csparams.sensmethod = sparams.sensmethod;
%       csparams.sensind = csdata.sensind;
%       csparams.sensmtx = csdata.sensmtx;
      
%     [x_zf, Tx_zf] = zero_filling(meas, csparams);
%     csparams.x0 = x_zf;

    % csparams.csmethod = 'PnP-ADMM'; % PnP-ADMM method
    %   csparams.rho = 13e2; % multiplier (noise regularier)
    %   csparams.denoiser = 'TV'; % total variation denoising
    %   csparams.maxiter = 100; % maximum number of iterations
    % 
    % x_tv = cs(sensmat, meas, csparams);
    % csparams.x0 = x_tv;
%  
    % csparams.csmethod = 'DAMP'; % D-AMP method (PnP-AMP)
    %   csparams.denoiser = 'TV'; % TV image denoising
    %   % csparams.denoiser = 'BM3D'; % BM3D image denoising
    %   % csparams.denoiser = 'FFDNet'; % FFDNet image denoising
    %   csparams.maxiter = 5; % maximum number of iterations
%     
%     csparams.channel_wise = false; % reconstruct each channel seperately
% %     csparams.csmethod = 'PnP-AP'; % PnP-AP method
% %     csparams.csmethod = 'PnP-IST'; % PnP-IST method
% %     csparams.csmethod = 'PnP-TwIST'; % PnP-TwIST method
% %     csparams.csmethod = 'PnP-FISTA'; % PnP-FISTA method
%     csparams.csmethod = 'PnP-ADMM'; % PnP-ADMM method
%       % csparams.rho = 1e-3; % multiplier (noise regularier)
%       % csparams.rho = 13e2; % multiplier (noise regularier)
%       csparams.rho = 2e2; % multiplier (noise regularier)
%     % csparams.csmethod = 'PnP-QCS'; % PnP-ADMM method
%     %   csparams.p = 3; % norm moment of BPDQ_p
%     %   csparams.rho = 12e2; % multiplier (noise regularier) [32x32x100%] (als) [0,12] lux
% %     csparams.csmethod = 'PnP-GAP'; % PnP-GAP method
% %     csparams.acc = false; % enable accelerated GAP
%       csparams.denoiser = 'TV'; % total variation denoising
%       csparams.maxiter = 50; % maximum number of iterations
% % 
%       % csparams.denoiser = 'BM3D'; % BM3D image denoising
%       % % csparams.denoiser = 'CBM3D'; % CBM3D color image denoising
%       % csparams.sigma   = [50 25 12  6]/MAXB; % noise deviation (to be estimated and adapted)
%       % csparams.maxiter = [10 10 10 10];
%       % % csparams.sigma   = [12  6]/MAXB; % noise deviation (to be estimated and adapted)
%       % % csparams.maxiter = [10 10];
% 
%       csparams.denoiser = 'FFDNet'; % FFDNet image denoising
%       if chnum == 1 || (isfield(csparams,'channel_wise') && csparams.channel_wise)  % grayscale image
%           load(fullfile('models','FFDNet_gray.mat'),'net');
%       else % color image
%           load(fullfile('models','FFDNet_color.mat'),'net');
%       end
%       csparams.net = vl_simplenn_tidy(net);
%       csparams.useGPU = true;
%       if csparams.useGPU
%           csparams.net = vl_simplenn_move(csparams.net, 'gpu') ;
%       end
%       csparams.ffdnetnorm_init = true; % use normalized video for the first 10 iterations
%       % csparams.ffdnetnorm = false; % normalize the video before FFDNet video denoising
%       csparams.ffdnetnorm = true; % normalize the video before FFDNet video denoising
%       csparams.sigma   = [50 25 12  6]/MAXB; % noise deviation (to be estimated and adapted)
%       csparams.maxiter = [10 10 10 10];
%       % csparams.sigma   = [12 6]/MAXB; % noise deviation (to be estimated and adapted)
%       % csparams.maxiter = [10 10];

% quant_delta = max(meas)/sparams.QUANT_STEPS;
% meas = meas/quant_delta;

t0 = tic;
sig_out = cs(sensmat,meas,csparams); % apply cs method
t_cs = toc(t0);

% sig_out = sig_out*quant_delta;
%
% sig_out = abs(sig_out);   % abs output
% sig_out = max(sig_out,0); % threshold output
samp_rc = reshape(sig_out,[rows,cols,chnum]);
recrmse = sqrt(immse(samp_rc,samp)); % root mean-square-error
recpsnr = psnr(samp_rc,samp);        % peak signal-to-noise-ratio
recssim = ssim(samp_rc,samp);        % structure similarity

% samp_nm = imnorm(samp_rc);
samp_rz = imresize(samp_rc,sparams.savesize,'nearest');
if isfield(csparams, 'denoiser')
    if isfield(sparams, 'transformshape')
        imwrite(samp_rz,sprintf('%s/%s%d_%s-%s_samprate%.2f_snr%ddb_%s-%s_rmse%.4f_psnr%.2f_ssim%.4f_t%.2f.png',outdir,...
            sparams.sampname,sparams.rows,sparams.sensmethod,sparams.transformshape,sparams.samprate,sparams.noisesnr,csparams.csmethod,csparams.denoiser,recrmse,recpsnr,recssim,t_cs));
    else
        imwrite(samp_rz,sprintf('%s/%s%d_%s_samprate%.2f_snr%ddb_%s-%s_rmse%.4f_psnr%.2f_ssim%.4f_t%.2f.png',outdir,...
            sparams.sampname,sparams.rows,sparams.sensmethod,sparams.samprate,sparams.noisesnr,csparams.csmethod,csparams.denoiser,recrmse,recpsnr,recssim,t_cs));
    end
else
    if isfield(sparams, 'transformshape')
        imwrite(samp_rz,sprintf('%s/%s%d_%s-%s_samprate%.2f_snr%ddb_%s_rmse%.4f_psnr%.2f_ssim%.4f_t%.2f.png',outdir,...
            sparams.sampname,sparams.rows,sparams.sensmethod,sparams.transformshape,sparams.samprate,sparams.noisesnr,csparams.csmethod,recrmse,recpsnr,recssim,t_cs));
    else
        imwrite(samp_rz,sprintf('%s/%s%d_%s_samprate%.2f_snr%ddb_%s_rmse%.4f_psnr%.2f_ssim%.4f_t%.2f.png',outdir,...
            sparams.sampname,sparams.rows,sparams.sensmethod,sparams.samprate,sparams.noisesnr,csparams.csmethod,recrmse,recpsnr,recssim,t_cs));
    end
end
figure; 
subplot(121);
  imshow(samp);
  title('Ground truth');
subplot(122);
  imshow(samp_rc);
  title(sprintf('%s recovery (PSNR %2.2fdB)',csparams.csmethod,recpsnr));

