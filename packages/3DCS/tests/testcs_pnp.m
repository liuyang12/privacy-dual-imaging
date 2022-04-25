%TESTCS test compressive sensing (CS) algorithms for sparse reconstruction
%   See also CS.
%% step 1. generate and save simulation dataset
clear; clc;
% close all;
% [1.0] parameters for simulation
sparams = []; % paramters for simulation
    sparams.rows = 64; % number of rows (m)
    sparams.cols = 64; % number of columns (n)
    sparams.samprate = 0.2; % sampling rate (gamma)
    sparams.noisesnr = 50; % signal-to-noise ratio (sigma)
    % sparams.noisesnr = inf; % signal-to-noise ratio (sigma)
    % sparams.sampname = 'gi'; % binary sample 
    % sparams.sampname = 'angrybird'; % binary sample 
    % sparams.sampname = 'phantom'; % simple gray-scale sample
    % sparams.sampname = 'cameraman'; % simple gray-scale sample
    sparams.sampname = 'peppers'; % simple gray-scale sample
    % sparams.sensmethod = 'binary'; % sensing matrix (multiple patterns) is binary
    % sparams.sensmethod = 'rademacher'; % sensing matrix (multiple patterns) is -1/+1
    sparams.sensmethod = 'mixhadamard'; % sensing matrix (multiple patterns) is binary
    sparams.SAMPLE_BINARY  = false; % sample is binary (1-Y, 0-N)
    sparams.savesize = 2; % saving size of the image
    % sparams.binsize = 1; % binning size of the spatial light modulator
REGENDATA = true; % flag of regenerating date (1-Y,0-N)
% REGENDATA = false; % flag of regenerating date (1-Y,0-N)
if REGENDATA
% [1] generate and save dataset
    gendata(sparams);
end

%% step 2. load and apply CS method to the dataset
% [2] load dataset
clear; clc;
% close all;
addpath('../utils')
addpath('../algorithms')
addpath(genpath('../packages'))

datadir    = '../../data';       % data directory
simdatadir = '../../data/sim/static';   % simulation data directory
outdir     = '../out';           % simulation output directory
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% [2.1] load simulation preferences
load(sprintf('%s/sim_prefs.mat',datadir));
% [2.2] load dataset
load(sprintf('%s/%s%d_samprate%.2f_snr%ddb.mat',simdatadir,...
    sparams.sampname,sparams.rows,sparams.samprate,sparams.noisesnr));
samp    = csdata.samp;    % sample
sensmat = csdata.sensmat; % sensing matrixs
meas    = csdata.meas;    % measurement vector

rows    = sparams.rows;
cols    = sparams.cols;
MAXB    = 255;           % maximum pixel value of the image (8-bit -> 255)

% [3] apply CS method
csparams = []; % parameters for CS method
    csparams.rows     = rows;
    csparams.cols     = cols;
    % csparams.srbasis  = 'haar'; % sparse representation basis
    % csparams.srbasis  = 'daubechies'; % sparse representation basis
    % csparams.srbasis  = 'dct'; % sparse representation basis

    % csparams.csmethod = 'gpsr'; % GPSR ell_1 solver
    % csparams.csmethod = 'ap'; % TV regularization method
    % csparams.csmethod = 'tv'; % TV regularization method
    % csparams.csmethod = 'gap'; % GAP ell_1 solver
    % csparams.csmethod = 'gap-tv'; % GAP-TV  solver
    % csparams.csmethod = 'gap2d'; % GAP2D solver

% [3.1] PnP-GAP-TV
    csparams.csmethod = 'pnp-admm'; % PnP-ADMM method
    csparams.rho = 1e-3; % multiplier (noise regularier)
%     csparams.csmethod = 'pnp-gaLiup'; % PnP-GAP method
%     csparams.acc = true; % enable accelerated GAP
      csparams.denoiser = 'tv'; % total variation denoising
      csparams.maxiter = 100; % maximum number of iterations
    sig_tv = pnp_admm(sensmat,meas,csparams);
   
% % [3.2] PnP-GAP-BM3D-TV    
%       csparams.x0 = sig_tv; % TV result as start point
%       csparams.denoiser = 'bm3d'; % BM3D image denoising
%       csparams.sigma   = [50 25 12  6]/MAXB; % noise deviation (to be estimated and adapted)
%       csparams.maxiter = [10 10 10 10];
    
% [3.3] PnP-GAP-WNNM-TV    
      csparams.x0 = sig_tv; % TV result as start point
      csparams.denoiser = 'wnnm'; % WNNM image denoising
      csparams.range   = 1; % signal value range during recovery
      csparams.sigma   = [12  6]/MAXB; % noise deviation (to be estimated and adapted)
      csparams.maxiter = [20 20];
    
% % [3.4] PnP-GAP-FFDNet-TV
%       csparams.x0 = sig_tv; % TV result as start point
%       csparams.denoiser = 'ffdnet'; % FFDNet image denoising
%       load(fullfile('models','FFDNet_gray.mat'),'net');
%       csparams.net = vl_simplenn_tidy(net);
%       csparams.useGPU = true;
%       if csparams.useGPU
%           csparams.net = vl_simplenn_move(csparams.net, 'gpu') ;
%       end
%       csparams.ffdnetvnorm_init = false; % use normalized video for the first 10 iterations
%       csparams.ffdnetvnorm = false; % normalize the video before FFDNet video denoising
%       % csparams.ffdnetvnorm = true; % normalize the video before FFDNet video denoising
%       csparams.sigma   = [25 12  6]/MAXB; % noise deviation (to be estimated and adapted)
%       csparams.maxiter = [10 10 10];

% sensmat = 2*sensmat-1;
% meas = meas - mean(meas);
tic;
sig_out = cs(sensmat,meas,csparams); % apply cs method
t_cs = toc;
% sig_out = abs(sig_out);   % abs output
% sig_out = max(sig_out,0); % threshold output
samp_rc = reshape(sig_out,rows,cols);
recrmse = sqrt(immse(samp_rc,samp)); % root mean-square-error
recpsnr = psnr(samp_rc,samp);        % peak signal-to-noise-ratio
recssim = ssim(samp_rc,samp);        % structure similarity

samp_nm = imnorm(samp_rc);
samp_rz = imresize(samp_nm,sparams.savesize,'nearest');
imwrite(samp_rz,sprintf('%s/%s%d_samprate%.2f_snr%ddb_%s_rmse%.4f_psnr%.2f_ssim%.4f_t%.1f.png',outdir,...
    sparams.sampname,sparams.rows,sparams.samprate,sparams.noisesnr,csparams.csmethod,recrmse,recpsnr,recssim,t_cs));
figure; 
imshow(samp_rz);


