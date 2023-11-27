%TESTCS3D Test compressive sensing (CS) algorithms for three dimensional
%sparse reconstruction.
%   See also CS3D.
%% step 1. generate and save simulation dataset
clear; clc;
% close all;
% [1.0] parameters for simulation
sparams = []; % paramters for simulation
    sparams.rows     = 32; % number of rows (m)
    sparams.cols     = 32; % number of columns (n)
    sparams.nframe   = 8;   % number of frames (F)
    sparams.samprate = 0.5; % sampling rate (gamma)
    sparams.noisesnr = 30;  % signal-to-noise ratio (sigma)
    sparams.sampname = 'foreman'; % grayscale video sample 
      sparams.width  = 176;   % width of the video
      sparams.height = 144;   % height of the video
      sparams.format = '420'; % YUV format ('420' for YUV 4:2:0 default)
    % sparams.sensmethod = 'binary'; % sensing matrix (multiple patterns) is binary
    sparams.sensmethod = 'mixhadamard'; % sensing matrix (multiple patterns) is binary
%     sparams.sensmethod = 'hadamard'; % sensing matrix (multiple patterns) is binary
%       sparams.transformshape = 'zigzag'; % sampling shape in transform domain
    sparams.SAMPLE_BINARY  = 0; % sample is binary (1-Y, 0-N)
    sparams.savesize = 5; % saving size of the image
    % sparams.binsize = 1; % binning size of the spatial light modulator
    sparams.SAME_SENSING_MATRIX = false; % same sensing matrix for each frame
% REGENDATA = true; % flag of regenerating date (1-Y,0-N)
REGENDATA = false; % flag of regenerating date (1-Y,0-N)
if REGENDATA
% [1] generate and save dataset
    gendata3d(sparams);
end

%% step 2. load and apply CS method to the dataset
% [2] load dataset
clear; clc;
% close all;
addpath('../utils')
addpath('../algorithms')
addpath(genpath('../packages'))

datadir    = '../../data';               % data directory
simdatadir = '../../data/sim/dynamic';   % simulation data directory
outdir     = '../vout';                  % simulation output directory
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% [2.1] load simulation preferences
load(sprintf('%s/sim_prefs.mat',datadir));
% [2.2] load dataset
load(sprintf('%s/%s%dby%d_samprate%.2f_snr%ddb.mat',simdatadir,...
    sparams.sampname,sparams.rows,sparams.nframe,sparams.samprate,...
    sparams.noisesnr));
samp    = csdata.samp;    % sample
sensmat = csdata.sensmat; % sensing matrixs
meas    = csdata.meas;    % measurement vector

rows    = sparams.rows;
cols    = sparams.cols;
nframe  = sparams.nframe;

% [3] apply CS method
cs3dparams = []; % parameters for CS method
    cs3dparams.rows     = rows;
    cs3dparams.cols     = cols;
    cs3dparams.nframe   = nframe;
%     cs3dparams.cs3dmethod = 'cs2d'; % 2D-CS for each frame
        % cs3dparams.srbasis  = 'haar'; % sparse representation basis
        % cs3dparams.srbasis  = 'daubechies'; % sparse representation basis
%         cs3dparams.srbasis  = 'dct'; % sparse representation basis

        % cs3dparams.csmethod = 'gpsr'; % GPSR ell_1 solver
        % cs3dparams.csmethod = 'tv'; % TV regularization method
%         cs3dparams.csmethod = 'gap'; % GAP ell_1 solver
%         cs3dparams.csmethod = 'gap-tv'; % GAP-TV  solver
        % cs3dparams.csmethod = 'gap2d'; % GAP2D solver
%     cs3dparams.cs3dmethod = 'tv3d'; % TV3D solver
    % cs3dparams.cs3dmethod = 'gap3d'; % GAP3D solver
    cs3dparams.cs3dmethod = 'gap_tv3d'; % GAP_TV3D solver
    cs3dparams.x0 = cs3d(sensmat,meas,cs3dparams); % apply cs method
    
%     cs3dparams.cs3dmethod = 'cs2d'; % CS2D solver
%         cs3dparams.csmethod = 'zero-filling'; % zero-filling for each subframe
%         cs3dparams.sensmethod = sparams.sensmethod;
%         cs3dparams.sensind = csdata.sensind;
%         cs3dparams.sensmtx = csdata.sensmtx;
%     cs3dparams.x0 = cs3d(sensmat,meas,cs3dparams); % apply cs method
%     cs3dparams.orthogonal_basis = true;
%     
    cs3dparams.cs3dmethod = 'pnp_admm3d'; % PnP_ADMM3D solver
        cs3dparams.rho = 1e-3;
        % cs3dparams.denoiser = 'vbm3d';
        % cs3dparams.denoiser = 'vbm4d';
        cs3dparams.denoiser = 'bm4d';
        cs3dparams.nosestim = false;
        cs3dparams.sigma = [20 10 5]/255;
        cs3dparams.maxiter = [5 5 5];

% profile clear
% profile -memory on
% tic;
% vout = cs3d(sensmat(:,:,1),meas,cs3dparams); % apply cs method
vout = cs3d(sensmat,meas,cs3dparams); % apply cs method
% t_cs = toc;
% profile report
% profile -memory off

%% [4] save all the frames
if strcmpi(cs3dparams.cs3dmethod,'cs2d')
    voutdir = sprintf('%s/%s%dby%d/%s',outdir,sparams.sampname,rows,nframe,cs3dparams.csmethod);
else
    voutdir = sprintf('%s/%s%dby%d/%s',outdir,sparams.sampname,rows,nframe,cs3dparams.cs3dmethod);
end
if isfield(cs3dparams, 'denoiser'), denoiser_suffix = ['-' cs3dparams.denoiser]; else denoiser_suffix = ''; end

if ~exist(voutdir,'dir')
    mkdir(voutdir);
end
for iframe = 1:nframe
    % sig_out = abs(sig_out);   % abs output
    % sig_out = max(sig_out,0); % threshold output
    samp_rc = reshape(vout(:,iframe),rows,cols);
    samp_nm = imnorm(samp_rc);
    samp_rz = imresize(samp_nm,sparams.savesize,'nearest');
    samp_or = reshape(samp(:,iframe),rows,cols); % original sample
    recrmse(iframe) = sqrt(immse(uint8(samp_nm*255),uint8(imnorm(samp_or)*255))); % root mean-square-error
    recpsnr(iframe) = psnr(uint8(samp_nm*255),uint8(imnorm(samp_or)*255));        % peak signal-to-noise-ratio
    recssim(iframe) = ssim(uint8(samp_nm*255),uint8(imnorm(samp_or)*255));        % structure similarity
%     imwrite(samp_rz,sprintf('%s/%s%d_%02d_samprate%.2f_snr%ddb_%s_rmse%.2f_psnr%.2f_ssim%.4f_t%.1f.png',voutdir,...
%         sparams.sampname,sparams.rows,iframe,sparams.samprate,sparams.noisesnr,cs3dparams.cs3dmethod,recrmse(iframe),recpsnr(iframe),recssim(iframe),t_cs));
%     imwrite(samp_rz,sprintf('%s/%s%d_%02d_samprate%.2f_snr%ddb_%s_rmse%.2f_psnr%.2f_ssim%.4f.png',voutdir,...
%         sparams.sampname,sparams.rows,iframe,sparams.samprate,sparams.noisesnr,cs3dparams.cs3dmethod,recrmse(iframe),recpsnr(iframe),recssim(iframe)));
end
figure; 
imshow(imresize(reshape(vout, rows, []),sparams.savesize));

imwrite(imresize(reshape(uint8(vout*255), rows, []),sparams.savesize),sprintf('%s/%s%d_samprate%.2f_snr%ddb_%s%s_rmse%.2f_psnr%.2f_ssim%.4f.png',voutdir,...
        sparams.sampname,sparams.rows,sparams.samprate,sparams.noisesnr,cs3dparams.cs3dmethod,denoiser_suffix,mean(recrmse),mean(recpsnr),mean(recssim)));

[mean(recrmse) mean(recpsnr) mean(recssim)]

% save the ground truth image sequence
imwrite(imresize(reshape(uint8(samp*255), rows, []),sparams.savesize),sprintf('%s/%s%d_groundtruth.png',voutdir,...
        sparams.sampname,sparams.rows));

