function [  ] = gendata3d( sparams )
%GENDATA3D generate and save three dimensioanal (3D) simulation dataset
%
% [0.0] add *utils*, *algorithms*, and *packages* to path 
addpath('../utils')
addpath('../algorithms')
addpath(genpath('../packages'))

% [0.1] directories configuration
datadir    = '../../data';                   % data directory
simdatadir = '../../data/sim/dynamic';       % simulation data directory
sampdir    = '../scenes/dynamic/noise_free'; % static sample directory
if ~exist(simdatadir,'dir')
    mkdir(simdatadir);
end

% [0.2] default paramter configuration
transformshape   = 'zigzag';     % sampled shape in transform domain
if isfield(sparams,'transformshape'), transformshape = sparams.transformshape; end

% [0.3] pre-calculation
switch lower(transformshape)
    case 'square' % center square in tramsform domain
        r = sqrt(sparams.samprate); 
    case 'zigzag' % center diamond in transform domain
        if sparams.samprate <= 1/2
            r = sqrt(2*sparams.samprate); 
        else
            r = 2 - sqrt(2*(1-sparams.samprate));
        end
    case 'circular' % center circle in transform domain
        if sparams.samprate <= pi/4
            r = 2*sqrt(1/pi*sparams.samprate); 
        else
            v = 0:1e-3:pi/2;
            x = (v+cos(v))./(1+sin(v));
            theta = interp1(x, v, sparams.samprate, 'linear', 'extrap'); % approximate
            r = 1/cos(pi/4-theta/2) + eps;
        end
    case 'mixed' % mixed orthogonal transform based sensing matrx
        %%% DO NOTHING HERE. %%%
    otherwise
        error('Unsupported shape %s in transform domain!',lower(sparams.transformshape));
end

% [1] load or generate data for simulation
nyqnum  = sparams.rows*sparams.cols; % number of Nyquist samping (N)
sampnum = round(sparams.samprate*nyqnum); % sampling number (M)
nframe  = sparams.nframe; % number of frames for simulation (F)
width   = sparams.width; % video width (W)
height  = sparams.height; % video height (H)
format  = sparams.format; % video format ('420' for YUV 4:2:0 default)
cols = sparams.cols;
rows = sparams.rows;

% [1.1] load video for simulation
samppath = sprintf('%s/%s.yuv',sampdir,sparams.sampname);
mov = yuv2mov(samppath,width,height,format); % read yuv format video
nmov = length(mov);
framespath = sprintf('%s/%s%d',sampdir,sparams.sampname,sparams.rows);
if ~exist(framespath,'dir')
    mkdir(framespath);
end
if nmov<nframe % number of frames in the video less than the requested 
               % number of frames
    error('not sufficient number of frames in %s.\n',samppath);
end
for imov = 1:nframe
    cursamp = mov(imov).cdata;
    if size(cursamp,3)>1
        cursamp = rgb2gray(cursamp);
    end
    if sparams.SAMPLE_BINARY % sample is binary
        level  = graythresh(cursamp);
        cursamp = im2bw(cursamp,level);
    end
    cursamp = imnorm(double(cursamp)); % normalization to fit grayscale sample
    if width>height % crop the middle square
        left  = floor((width-height)/2)+1;
        right = floor((width+height)/2);
        cursamp_sq = cursamp(:,left:right);
    else
        left  = floor((height-width)/2)+1;
        right = floor((height+width)/2);
        cursamp_sq = cursamp(left:right,:);
    end
    cursamp_sq = imresize(cursamp_sq,[sparams.rows sparams.cols]);
    samp(:,imov) = cursamp_sq(:); % all the samples in the video (m*n)*F
    imwrite(imresize(cursamp_sq,sparams.savesize,'nearest'), sprintf('%s/frame%04d.png', framespath, imov));
end

if sparams.SAME_SENSING_MATRIX
    sensframe = 1;
else
    sensframe = nframe;
end
% [1.2] generate sensing matrix and get measurement vector as dataset
switch sparams.sensmethod
    case 'binary' % binary pseudo-random sensing matrix (A)
        sensmat = rand(sampnum,nyqnum,sensframe)>0.5; % binary sensing matrix (A)
    case 'gaussian' % Gaussian sensing matrix (A)
        sensmat = randn(sampnum,nyqnum,sensframe); 
        % sensmat = 1/sampnum*randn(sampnum,nyqnum); % Gaussian sensing matrix (A)
    case {'fourier', ... % Discrete Fourier transform (DFT) sensing matrix [-1/i,1/i] [orthogonal] 'fourier'='dft'
          'dft', ...     % Discrete Fourier transform (DFT) sensing matrix [-1/i,1/i] [orthogonal] 'fourier'='dft'
          'dct', ...     % Discrete cosine transform (DCT) sensing matrix [-1,1] [orthogonal] 'dct'='cosine'
          'cosine', ...  % Discrete cosine transform (DCT) sensing matrix [-1,1] [orthogonal] 'dct'='cosine'
          'haar', ...    % Haar wavelet transform sensing matrix [-1,1] [orthogonal] 
          'walsh', ...   % Walsh-Hadamard sensing method -1/+1 [orthogonal] 'hadamard'='walsh'
          'hadamard'}    % Walsh-Hadamard sensing method -1/+1 [orthogonal] 'hadamard'='walsh'
        sparams.SAME_SENSING_MATRIX = true;
        sparams.orthogonal_basis = true;
        x = 1:cols;
        y = 1:rows;
        [X,Y] = meshgrid(x,y);
        rad = eye(nyqnum);
        if sum(strcmpi(sparams.sensmethod,{'fourier','dft'})) % Fourier zero-frequency at the center
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
        switch lower(sparams.sensmethod)
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
    case 'mixhadamard' % mixed Hadamard sensing method
        Hr = hadamard(sparams.rows); % Hadamard matrix in rows
        Hc = hadamard(sparams.cols); % Hadamard matrix in columns
        Had = kron(Hc',Hr');
        for iframe = 1:sensframe
            Rad = double(rand(nyqnum,1)>0.5)*2-1; % Rademacher distribution
            allrowset = randperm(nyqnum);
            rowset = allrowset(randperm(sampnum));
            sensmat(:,:,iframe) = Had(rowset,:)*diag(Rad);
            csdata.sensmtx = Had;
            csdata.sensind = rowset;
            csdata.Rad     = Rad;
        end
    otherwise
        error('unsupported sensing method.\n');
end
sensmat = squeeze(sensmat);
meas_nonoise = fm_3dcs(sensmat,samp); % forward model of 3D-CS
noisenorm  = randn(sampnum,nframe);
noise = norm(meas_nonoise,2)/norm(noisenorm,2)*10^(-sparams.noisesnr/20)*noisenorm;
meas = meas_nonoise+noise; % measurements (y)

% [1.3] save generated dataset for simulation
csdata.samp    = samp;
csdata.sensmat = sensmat;
csdata.meas    = meas;
save(sprintf('%s/%s%dby%d_samprate%.2f_snr%ddb.mat',simdatadir,...
    sparams.sampname,sparams.rows,sparams.nframe,sparams.samprate,...
    sparams.noisesnr),'csdata');
% [1.4] save sparams as simulation preferences
save(sprintf('%s/sim_prefs.mat',datadir),'sparams');

end

