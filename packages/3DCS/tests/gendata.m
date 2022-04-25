function [  ] = gendata( sparams )
%GENDATA generate and save simulation dataset
%
% [0.0] add *utils*, *algorithms*, and *packages* to path 
addpath('../utils')
addpath('../algorithms')
addpath(genpath('../packages'))

% [0.1] directories configuration
datadir    = '../../data';       % data directory
simdatadir = '../../data/sim/static';   % simulation data directory
sampdir    = '../scenes/static'; % static sample directory
if ~exist(simdatadir,'dir')
    mkdir(simdatadir);
end

% [0.2] default paramter configuration
transformshape   = 'zigzag';     % sampled shape in transform domain
SYMMETRY         = true;         % Hermition (conjugate) symmetry of (discrete) Fourier transform 
NUMPHSHIFT       = 4;            % number of phase shifts to get the complex value of the DFT sensing matrix (3 as alternative)
QUANTIZATION     = false;        % apply quantization to the direct measurement
  QUANT_STEPS    = 255;          % steps of quantization
dither_num       = 1;            % number of dithered data points per measurement

if isfield(sparams,'transformshape'), transformshape = sparams.transformshape; end
if isfield(sparams,'SYMMETRY'),             SYMMETRY = sparams.SYMMETRY;       end
if isfield(sparams,'NUMPHSHIFT'),         NUMPHSHIFT = sparams.NUMPHSHIFT;     end
if isfield(sparams,'QUANTIZATION'),     QUANTIZATION = sparams.QUANTIZATION;   end
if isfield(sparams,'QUANT_STEPS'),       QUANT_STEPS = sparams.QUANT_STEPS;    end
if isfield(sparams,'dither_num'),         dither_num = sparams.dither_num;    end

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
rows = sparams.rows;
cols = sparams.cols;
nyqnum  = rows*cols; % number of Nyquist samping (N)
sampnum = round(sparams.samprate*nyqnum); % sampling number (M)

% [1.1] load sample for simulation
samppath = sprintf('%s/%s',sampdir,sparams.sampname);
samp = imreadallfmt(samppath); % read sample image
chnum = size(samp, 3); % number of color channels (1 grayscale, 3 for color)
if chnum>1 && (~isfield(sparams,'SAMPLE_RGB') || ~sparams.SAMPLE_RGB)
    samp = rgb2gray(samp);
    chnum = 1;
end
if sparams.SAMPLE_BINARY % sample is binary
    level = graythresh(samp);
    samp = im2bw(samp, level);
end
% samp = imnorm(samp); % normalization to fit grayscale sample
samp = double(samp)/max(double(samp(:))); % normalized to [0,1]
if isfield(sparams,'SIMREALSIZE') && sparams.SIMREALSIZE
    samp_rz = imresize(samp, sparams.cropsize);
    samp_rzvec = reshape(samp_rz, [sparams.cropsize(1)*sparams.cropsize(2) chnum]);
end
samp = imresize(samp, [rows cols]);

imwrite(imresize(samp,sparams.savesize,'nearest'), sprintf('%s/%s%d.png', sampdir, sparams.sampname, rows));

sparams.orthogonal_basis = false; % [indicator] orthogonal basis or not
% [1.2] generate sensing matrix and get measurement vector as dataset
switch lower(sparams.sensmethod)
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
        error('unsupported sensing method %s.\n', sparams.sensmethod);
end

if QUANTIZATION
    if sum(strcmpi(sparams.sensmethod,{'rademacher','dct','cosine','haar','walsh','hadamard','mixhadamard'})) % [-1,1]
        pmeas_nonoise = (1+sensmat)/2*reshape(samp, [nyqnum chnum]); % (1+H)/2 [0,1]
        nmeas_nonoise = (1-sensmat)/2*reshape(samp, [nyqnum chnum]); % (1-H)/2 [0,1]
        
        pnmeas_nonoise = [pmeas_nonoise, nmeas_nonoise];
        
        % add noise according to measurement SNR
        noisenorm  = randn(sampnum,chnum*2,dither_num);
        noise = norm(pnmeas_nonoise(:),2)/norm(noisenorm(:),2)*10^(-sparams.noisesnr/20)*noisenorm;
        pnmeas = pnmeas_nonoise+noise; % measurements (y)
        
        % quantization after adding measurement noise
        % VMAX = max(pnmeas(:)); % use the measurement maximum
        VMAX = max(sum(reshape(samp, [nyqnum chnum]),1)); % overall measurement maximum 
        pnmeas = round(pnmeas/VMAX*QUANT_STEPS)/QUANT_STEPS*VMAX;
        
        % average over all temporally dithered data points
        pnmeas = mean(pnmeas,3);
        
        pmeas = pnmeas(:,1:chnum,:);
        nmeas = pnmeas(:,chnum+1:end,:);
        
        % final measurements reflecting the sensing matrix
        meas = pmeas-nmeas;
        
    elseif sum(strcmpi(sparams.sensmethod,{'fourier','dft'})) % (Discrete) Fourier transform [-1/-i,1/i]
        % TODO: fix for width and height being even numbers when sampling
        % rate over pi/4.
        if SYMMETRY % empploying Hermitian (conjugate) symmetry
            % C = ceil((sampnum+1)/2); % center in 1D
            c0 = (x0-1)*rows + y0; % center in 1D at the original transform domain
            C = find(sensind==c0); % center in 1D at the sampled transform domain
        else
            C = sampnum;
        end
        
        % phase shifting to get the complex value of the Fourier sensing
        % matrix
        phmeas_nonoise = zeros([C, chnum, NUMPHSHIFT]);
        for iph = 1:NUMPHSHIFT
            phmeas_nonoise(:,:,iph) = 1/2*(1+cos(angle(sensmat(1:C,:))+(iph-1)*2*pi/NUMPHSHIFT))*reshape(samp, [nyqnum chnum]); % I_phi
        end
        
        % add noise according to measurement SNR
        noisenorm = randn([C, chnum, NUMPHSHIFT, dither_num]);
        noise = norm(phmeas_nonoise(:),2)/norm(noisenorm(:),2)*10^(-sparams.noisesnr/20)*noisenorm;
        phmeas = phmeas_nonoise+noise; % measurements (y)
        
        % quantization after adding noise
        % VMAX = max(phmeas(:)); % measurement maximum 
        VMAX = max(sum(reshape(samp, [nyqnum chnum]),1)); % overall measurement maximum 
        phmeas = round(phmeas/VMAX*QUANT_STEPS)/QUANT_STEPS*VMAX;
        
        % average over all temporally dithered data points
        phmeas = mean(phmeas,4);
        
        %  measurements reflecting the sensing matrix
        switch NUMPHSHIFT
            case 4 % four-step phase shifting [0,pi/2,pi,3*pi/2] -> working
                % e^(j*a) = cos(a)+j*sin(a) = (I0-I2) + j*(I3-I1)
                meas =      (phmeas(:,:,1)-phmeas(:,:,3)) ...
                       + 1i*(phmeas(:,:,4)-phmeas(:,:,2)); 
            case 3 % three-step phase shifting [0,2*pi/3,4*pi/3] -> working
                % e^(j*a) = cos(a)+j*sin(a) = 2/3*(2*I0-I1-I2)+j*2^√3/3*(I2-I1)
                meas = 2/3*(2*phmeas(:,:,1)-phmeas(:,:,2)-phmeas(:,:,3)) ...
                       + 1i*2*sqrt(3)/3*(phmeas(:,:,3)-phmeas(:,:,2));
            case 2 % two-step phase shifting [0,pi/2] -> not working well
                % e^(j*a) = cos(a)+j*sin(a) = 2*(I0-I_)-j*2*(I1-I_)
                Imean = 1/2*ones([C nyqnum])*reshape(samp, [nyqnum chnum]);
                meas = 2*(phmeas(:,:,1)-Imean) - 1i*2*(phmeas(:,:,2)-Imean);
            otherwise
                error('Unsupported number of phase shifts %d in %s sensing matrix',NUMPHSHIFT,sparams.sensmethod);
        end

        if SYMMETRY % empploying Hermitian (conjugate) symmetry
            if r<1/2 || (mod(rows,2)&&mod(cols,2)) % sampled space symmetric at the center frequency
                meas = [meas;conj(meas(end-1:-1:1,:))];
            else % sampled space not symmetric at the center frequency
                TD_vec = zeros([rows*x0 chnum]);
                TD_vec(sensind(1:C),:) = meas;
                TD_vec(sensind(C+(1:rows-y0),:)) = conj(meas(end-1:-1:end-(rows-y0),:));
                TD = reshape(TD_vec, [rows x0, chnum]);
                % TD = reshape([meas;conj(meas(end-1:-1:end-(rows-y0),:))], [rows,x0,chnum]);
                if mod(rows,2) % odd
                    TDsym = conj(TD(end:-1:1, end-1:-1:2-mod(cols,2),:));
                else % even
                    TDsym = conj(TD([1 end:-1:2], end-1:-1:2-mod(cols,2),:));
                end
                TD = cat(2,TD,TDsym);
                TD_vec = reshape(TD, [nyqnum chnum]);
                meas = TD_vec(sensind,:);
            end
        end
        
    elseif sum(strcmpi(sparams.sensmethod,{'binary','video'})) % [0,1]
        if ~ (isfield(sparams,'SIMREALSIZE') && sparams.SIMREALSIZE)
            meas_nonoise = sensmat*reshape(samp, [nyqnum chnum]);
        end
        
        % add noise according to measurement SNR
        noisenorm  = randn(sampnum,chnum,dither_num);
        noise = norm(meas_nonoise,2)/norm(noisenorm,2)*10^(-sparams.noisesnr/20)*noisenorm;
        meas_unquant = meas_nonoise+noise; % unquantized measurements (y)
        
        % quantization after adding measurement noise
        VMAX = max(meas_unquant(:)); % measurement maximum
        % VMAX = max(sum(reshape(samp, [nyqnum chnum]),1)); % overall measurement maximum 
        meas = round(meas_unquant/VMAX*QUANT_STEPS)/QUANT_STEPS*VMAX;
        
        % average over all temporally dithered data points
        meas = mean(meas,3);
        
    else % Gaussian sensing matrix
        
    end
    
    VMAX
else
    if ~ (isfield(sparams,'SIMREALSIZE') && sparams.SIMREALSIZE)
        meas_nonoise = sensmat*reshape(samp, [nyqnum chnum]);
    end
    noisenorm  = randn(sampnum,chnum,dither_num);
    noise = norm(meas_nonoise,2)/norm(noisenorm,2)*10^(-sparams.noisesnr/20)*noisenorm;
    meas = meas_nonoise+noise; % measurements (y)
    
    % average over all temporally dithered data points
    meas = mean(meas,3);
end

% [1.3] save generated dataset for simulation
csdata.samp    = samp;
csdata.sensmat = sensmat;
csdata.meas    = meas;
save(sprintf('%s/%s%d_%s_samprate%.2f_snr%ddb.mat',simdatadir,sparams.sampname,rows,sparams.sensmethod,sparams.samprate,sparams.noisesnr),'csdata','-v7.3');
% [1.4] save sparams as simulation preferences
save(sprintf('%s/sim_prefs.mat',datadir),'sparams');

end

