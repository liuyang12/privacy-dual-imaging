%FIG02_RECON_PROCEDURE_TOUCH 
% Figure 2 | Inversion procedure for revealing imaging privacy threats from
% an ambient light sensor.
clear; clc;
close all;
% [0] environment configuration
addpath(genpath('../packages')); % packages
addpath('../utils'); % utilities

rawdata_dir = '../rawdata/'; % raw data directory
dataset_dir = '../dataset/'; % dataset directory
results_dir = '../results/'; % results

scfile      = 2011; % scene file(s), average if multiple
bgfile      = []; % background file(s), average if multiple

measidx = [28]; % index for plotting interested frames
% rhos = [3]*1e2;
ges_suffix = '_touch_fig';

GIFFLAG = false;
  gif_delaytime = 0.5; % [s] 
csmethod = 'pnp';
  denoiser = 'ffdnet';
  
nfigcol = 1;

% [1] parameters
SAMPLE    = 'touch-hand'; % sample (scene) name
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
HFLIP     = false; % horizontal flip
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
params.interfactor = 1.5; % interval factor (to determine the segmentation)
params.measthresh = 'mixed'; % single channel [ambient light sensor output]
params.processmeth = 'mean'; % processing method to obtain each measurement [maximum of a period]
params.avginterval = [0.2, 0.8];
% params.processmeth = 'max'; % prtocessing method to obtain each measurement [maximum of a period]
% params.avginterval = [0, 1];
% params.avginterval = [-0.2, 1.2];
%   params.patseg = patseg;

% magsize = 16; % magnification factor
magsize = [256 256]; % magnified size

if ~exist([dataset_dir SAMPLE], 'dir')
    mkdir([dataset_dir SAMPLE]);
end
if ~exist([results_dir SAMPLE],'dir')
    mkdir([results_dir SAMPLE]);
end

ifile = scfile;

vname = sprintf('%s%dx%d_%dcm_%dfps_%03d',SAMPLE,h,v,distance,framerate,ifile);
datafile = [dataset_dir SAMPLE '/' vname '.mat'];

if exist(datafile, 'file') && (exist('RERAWREAD','var') && ~RERAWREAD)
    load(datafile,'measall');
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
        % table = readtable(csvfname,'Format','%f%f%d%s','VariableNamingRule','preserve');
        table = readtable(csvfname,'Format','%f%d%s','VariableNamingRule','preserve');
        Iall = table2array(table(:,1));
        % % RGBW raw data
        % table = readtable(csvfname,'Format','%f%f%f%f%f%s','VariableNamingRule','preserve');
        % Iall = table2array(table(:,1:3))';
    else
        error('raw data %s not found.',vname);
    end
    dlen = length(Iall)/nColor;
    I=reshape(Iall,nColor,[])';

    % I = round((I-28.96)/1.18);

    % % [3] turn the raw measurements to the corresponding correlated
    % measurements
    [measrawall patsegall] = readvidmeasdata(I,params); 

    if iscell(measrawall)
        if DIFFMEAS % differential measurements
            for i = 1:length(measrawall)
                measall{i} = measrawall{i}(:,1:2:end) - measrawall{i}(:,2:2:end);
            end
        else
            measall = measrawall;
        end
        
        if COLORPROJ % color pattern projection
            for i = 1:length(measall)
                measall{i} = reshape(measall{i},[],nch)';
            end
        end
    else
        if DIFFMEAS % differential measurements
            measall = measrawall(:,1:2:end) - measrawall(:,2:2:end);
        else
            measall = measrawall;
        end
        
        if COLORPROJ % color pattern projection
            measall = reshape(measall,[],nch)';
        end
        measall = {measall};
    end

    % save raw data
    save(datafile,'measall');
end

if ~exist('measidx', 'var') || isempty(measidx)
    measidx = 1:length(measall);
end
nfigrow = ceil(length(measidx) / nfigcol);
if FIGFLAG, fig = figure; end

h_separator = 5; % Width of the horizontal separator in pixels
v_separator = 5; % Height of the vertical separator in pixels

stacked_image_width = nfigcol * magsize(2) + (nfigcol - 1) * h_separator;
stacked_image_height = nfigrow * magsize(1) + (nfigrow - 1) * v_separator;

stacked_image = ones(stacked_image_height, stacked_image_width, nColor);

nframe = length(measidx);
meas3d = reshape(cell2mat(measall(measidx)), [], nframe); % [M,F]

% for i = 1:length(measidx)
%     meas = measall{measidx(i)};

%% [4] reconstruction of the scene, color channel by color channel
%       the denoising-based approximate message passing (D-AMP) algorithm
%       is employed here
% [4.1] load the sensing matrix 
MAXB = 255;
chnum = 1;

para.rows = v; % [vertical] height
para.cols = h; % [horizontal] width
para.nframe = nframe; % number of frames

samprate = round(size(meas3d,1)/(h*v)/0.25)*0.25; % sampling rate 

load(sprintf('../packages/data/sim/static/peppers%d_hadamard_samprate%.2f_snr50db.mat',h,samprate),'csdata');
% load(sprintf('../packages/data/sim/static/mit_logo%d_video_samprate%.2f_snr50db.mat',h,samprate),'csdata');

% load(sprintf('%s/%s%d_samprate%.2f_snr%ddb.mat',simdatadir,...
%     sparams.sampname,sparams.rows,sparams.samprate,sparams.noisesnr));

% reconstruction

para.cs3dmethod = 'cs2d'; % 2D-CS for each frame

if ~exist('csmethod', 'var')
    csmethod = 'zero-filling';
end
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
        % x_zf = cs3d(csdata.sensmat, meas3d, para);
        [x_zf, Tx_zf] = zero_filling(meas3d, para);
        iminit = reshape(x_zf,[v h]);
        
        para.x0v = x_zf; % zero-filling as initialization
        
        % para.csmethod = 'DAMP'; % D-AMP method (PnP-AMP)
        %   para.denoiser = 'FFDNet'; % FFDNet image denoising
        %   % para.denoiser = 'BM3D'; % FFDNet image denoising
        %   para.maxiter = 5; % maximum number of iterations
        
        para.orthogonal_basis = true;
        
        para.channel_wise = true; % channel-wise reconstruction
        % para.csmethod = 'PnP-ADMM'; % PnP-ADMM method
        %   para.rho = 3e2; % multiplier (noise regularier) [32x32x100%] (als) [0,18] lux
        para.csmethod = 'PnP-QCS'; % PnP-QCS method
          para.p = 3; % norm moment of BPDQ_p
          para.rho = 6e2; % multiplier (noise regularier) [32x32x100%] (als) [0,18] lux
          if exist('rhos', 'var')
              para.rhos = rhos;
          end
        % para.csmethod = 'PnP-QCS'; % PnP-ADMM method
        %   para.p = 2; % norm moment of BPDQ_p
        % para.rho = 0.1e2; % multiplier (noise regularier) [32x32x100%] (als) [0,18] lux
        
        % para.csmethod = 'PnP-GAP'; % PnP-GAP method
        % para.acc = false; % enable accelerated GAP
        
        switch lower(denoiser)
            case 'tv'
                para.denoiser = 'TV'; % total variation denoising
                para.maxiter =  50; % maximum number of iterations
                para.maxiter = 100; % maximum number of iterations
            case 'bm3d'
                para.denoiser = 'BM3D'; % BM3D image denoising
                % para.denoiser = 'CBM3D'; % color BM3D image denoising
                % para.sigma   = [12  6]/MAXB; % noise deviation (to be estimated and adapted)
                % para.maxiter = [10 10];
                para.sigma   = [ 6]/MAXB; % noise deviation (to be estimated and adapted)
                para.maxiter = [10];
            case 'wnnm'
                para.denoiser = 'WNNM'; % WNNM image denoising
                para.range   = 1; % signal value range during recovery
                para.sigma   = [12  6]/MAXB; % noise deviation (to be estimated and adapted)
                para.maxiter = [20 20];
            case 'ffdnet'
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
                para.ffdnetnorm = true; % normalize the video before FFDNet video denoising
                para.sigma   = [12  6]/MAXB; % noise deviation (to be estimated and adapted)
                para.maxiter = [10 10];
            otherwise
                disp('unsupported denoiser.');
        end
    case 'pnp3d'
        para.cs3dmethod = 'cs2d'; % CS2D solver
            para.csmethod = 'zero-filling'; % zero-filling for each subframe
            para.sensmethod = 'hadamard';
            para.sensind = csdata.sensind;
            para.sensmtx = csdata.sensmtx;
        para.x0 = cs3d(csdata.sensmat,meas3d,para); % apply cs method
        para.orthogonal_basis = true;

        para.cs3dmethod = 'pnp_admm3d'; % PnP_ADMM3D solver
            para.rho = 5e2;
            cs3dparams.denoiser = denoiser;
            para.nosestim = false;
            para.sigma = [100 50 20 10]/255;
            para.maxiter = [5 5 5 5];
    otherwise
        disp('unsupported csmethod.');
end

tstart = tic;
x3d = cs3d(csdata.sensmat, meas3d, para); % apply cs method [integrated channels]
t_cs = toc(tstart);

for i = 1:nframe
    x = x3d(:,i);
% if COLORPROJ % color projection
%     x = x*inv(colormtx)';
% end

imraw = reshape(x, [v h chnum]);

% apply intensity compensation to the recovered image
if BRCORRFLAG
    a = 2;
    % b = 0.7;
    b = 0.8;
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
    iminit = iminit.*(1+((X-x0).^2/a^2+(Y-y0).^2/b^2)/d^2);
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

imrecon = imnorm(imraw);
iminit = imnorm(iminit);
if HFLIP % flip horizontally
    imrecon = flip(imrecon,2);
    iminit = flip(iminit,2);
end
if UPSIDEDOWN % flip horizontally
    imrecon = rot90(imrecon,2);
    iminit = rot90(iminit,2);
end
if nColor > 1 || (COLORPROJ && nch > 1)
    blank = zeros(size(imrecon,[1 2]));
    red_channel = cat(3, imrecon(:,:,1), blank, blank);
    green_channel = cat(3, blank, imrecon(:,:,2), blank);
    blue_channel = cat(3, blank, blank, imrecon(:,:,3));

    if FIGFLAG
        fig = figure;
        subplot(221); imshow(imresize(red_channel,magsize,'nearest')); title('Red channel');
        subplot(222); imshow(imresize(green_channel,magsize,'nearest')); title('Green channel');
        subplot(223); imshow(imresize(blue_channel,magsize,'nearest')); title('Blue channel');
        subplot(224); imshow(imresize(imrecon,magsize,'nearest')); title('RGB color');
    end
else
    if FIGFLAG
        % fig = figure;
        subplot(nfigrow, nfigcol, i);
        imshow(imresize(imrecon,magsize,'nearest')); title('Reconstructued');
        title(sprintf('Frame #%d',i));
    end
end

% save each frame of the reconstructed sequence to the stacked image
row_idx = floor((i-1) / nfigcol);
col_idx = mod(i-1, nfigcol);

start_row = row_idx * (magsize(1) + v_separator) + 1;
end_row = start_row + magsize(1) - 1;
start_col = col_idx * (magsize(2) + h_separator) + 1;
end_col = start_col + magsize(2) - 1;

stacked_image(start_row:end_row, start_col:end_col, :) = imresize(imrecon,magsize,'nearest');

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

if ~exist([results_dir, SAMPLE, '/', vname], 'dir')
    mkdir([results_dir, SAMPLE, '/', vname]);
end

if nColor > 1 || (COLORPROJ && nch > 1), color_suffix = '_color'; else color_suffix = ''; end
if isfield(para, 'denoiser'), denoiser_suffix = ['-' para.denoiser]; else denoiser_suffix = ''; end
if BRCORRFLAG, brcorr_suffix = '_brcorr'; else brcorr_suffix = ''; end
if ~exist('ges_suffix', 'var'), ges_suffix = ''; end

savename = sprintf('%s%s/%s/%s_%s%s%s%s%s', results_dir, SAMPLE, vname, vname, para.csmethod, denoiser_suffix, color_suffix, brcorr_suffix, ges_suffix);

% if nColor > 1 || (COLORPROJ && nch > 1)
%     if isfield(para, 'denoiser')
%         if BRCORRFLAG
%             savename = sprintf('%s%s/%s/%s_%s-%s_color_brcorr', results_dir, SAMPLE, vname, vname, para.csmethod, para.denoiser);
%         else
%             savename = sprintf('%s%s/%s/%s_%s-%s_color', results_dir, SAMPLE, vname, vname, para.csmethod, para.denoiser);
%         end
%     else
%         if BRCORRFLAG
%             savename = sprintf('%s%s/%s/%s_%s_color_brcorr', results_dir, SAMPLE, vname, vname, para.csmethod);
%         else
%             savename = sprintf('%s%s/%s/%s_%s_color', results_dir, SAMPLE, vname, vname, para.csmethod);
%         end
%     end
% else
%     if isfield(para, 'denoiser')
%         if BRCORRFLAG
%             savename = sprintf('%s%s/%s/%s_%s-%s_brcorr', results_dir, SAMPLE, vname, vname, para.csmethod, para.denoiser);
%         else
%             savename = sprintf('%s%s/%s/%s_%s-%s', results_dir, SAMPLE, vname, vname, para.csmethod, para.denoiser);
%         end
%     else
%         if BRCORRFLAG
%             savename = sprintf('%s%s/%s/%s_%s_brcorr', results_dir, SAMPLE, vname, vname, para.csmethod);
%         else
%             savename = sprintf('%s%s/%s/%s_%s', results_dir, SAMPLE, vname, vname, para.csmethod);
%         end
%     end
% end

savename_stackedimage = sprintf('%s%s/%s_%s%s%s%s%s', results_dir, SAMPLE, vname, para.csmethod, denoiser_suffix, color_suffix, brcorr_suffix, ges_suffix);
if exist('GIFFLAG', 'var') && GIFFLAG
    gifname = [savename_stackedimage '.gif'];
    imfinal = uint8(imresize(imrecon,magsize,'nearest')*255);
    if i == 1
        imwrite(imfinal, gifname, 'gif', 'LoopCount', Inf, 'DelayTime', gif_delaytime);
    else
        imwrite(imfinal, gifname, 'gif', 'WriteMode', 'append', 'DelayTime', gif_delaytime);
    end
end

savename = sprintf('%s_frame%02d', savename, measidx(i));

imwrite(imresize(imrecon,magsize,'nearest'),[savename '.png']);
end

figure; imshow(stacked_image);
imwrite(stacked_image,[savename_stackedimage '_stacked.png']);

if FIGFLAG, saveas(fig,[savename '.fig']); end


%% plot the raw measurements and results for demonstration of the reconstruction procudure

figdir = './savedfig_touch';
if ~exist(figdir, 'dir')
    mkdir(figdir);
end

% change default font name and font size
textfontsize = 20; % 14
legendfontsize = 20; % 14
linewidth = 4; % 2

% Change default axes fonts.
set(0,'DefaultAxesFontName','Myriad Pro');
% set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultAxesFontSize',textfontsize);
set(0,'DefaultAxesFontWeight','normal');
% Change default text fonts.
set(0,'DefaultTextFontname','Myriad Pro');
% set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontSize',textfontsize);
set(0,'DefaultTextFontWeight','normal');

meas = meas3d;

% Figure 1(a) | raw measurements
fig1 = figure('position',[50,100,425,250]);
set(gca, 'Color', 'white');
plot(meas,'b-','linewidth',linewidth);
xlabel('Sampling index over time'); xlim([0,length(meas)+1]);
ylabel('Measurement intensity'); ylim([-5,25]); yticks(-5:10:25);
set(gca,'linewidth',linewidth);

saveas(fig1, [figdir '/fig01_a_recon_procedure_rawmeas.svg']);

% Figure 1(b) | transform domain
para.csmethod = 'zero-filling'; % naive zero-filling (not CS though)
                                %  works only with orthogonal bases
para.sensmethod = 'hadamard';
para.sensind = csdata.sensind;
para.sensmtx = csdata.sensmtx;
[x_zf, Tx_zf] = zero_filling(meas, para);

fig2 = figure('position',[50,100,300,250]);
set(gca, 'Color', 'white');
imagesc(reshape(Tx_zf,[v h])); axis image; 
  c = colorbar; caxis([-5 25]); c.Ticks=-5:10:25; c.LineWidth=linewidth;
  xticks(10:10:30);yticks(10:10:30);
set(gca,'linewidth',linewidth);

saveas(fig2, [figdir '/fig01_b_recon_procedure_transform_domain.svg']);

% Figure 1(c) | raw initialization
fig3 = figure('position',[50,100,300,250]);
set(gca, 'Color', 'white');
imagesc(iminit); axis image; colormap('gray');
  c = colorbar; caxis([0 1]); c.Ticks=0:1; c.LineWidth=linewidth;
  xticks(10:10:30);yticks(10:10:30);
set(gca,'linewidth',linewidth);

saveas(fig3, [figdir '/fig01_c_recon_procedure_raw_initialization.svg']);

% Figure 1(d) | final reconstruction
fig4 = figure('position',[50,100,300,250]);
set(gca, 'Color', 'white');
imagesc(imrecon); axis image; colormap('gray');
  c = colorbar; caxis([0 1]); c.Ticks=0:1; c.LineWidth=linewidth;
  xticks(10:10:30);yticks(10:10:30);
set(gca,'linewidth',linewidth);

saveas(fig4, [figdir '/fig01_d_recon_procedure_raw_final_recon.svg']);

