function [ sig_out ] = cs( sensmat,meas,params )
%CS compressive sensing algorithms in solving sparse reconstruction problem
%   out=CS(sensmat,meas,params) returns the N-by-1 vector of the sparse 
%   reconstruction result, where sensmat is the M-by-N sensing matrix or 
%   measurement matrix, M is the number of measurements and N is the 
%   Nyquist number of the orignal signal, meas is the M-by-1 measurement 
%   vector that satisies the compressiver sampling equation 
%   meas=sensmat*sig_out+nois, and nois here is the measurement noise on 
%   meas, params is a structure specifying the parameters of the solver, 
%   including the algorithm (method) with the corresponding options.
%   See also GPSR, AP, TV, TVAL3, GAP_TV.
% 
%   Copyright(C) 2016-2017 by <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>.
%   Last modified Feb 9, 2017.

% % [0] add *algorithms* and *packages* to operation path, might be done out
% % of this function block
% addpath('../algorithms')
% addpath(genpath('../packages'))

% [1] basic parameters
opts = []; % options for various algorithms

% differential sensing matrix (row subtraction) and the corresponding
% measurements
if isfield(params, 'DIFFSENS') && params.DIFFSENS
    sensmat = diff(sensmat, 1);
    meas    = diff(meas, 1);
end

% [2] apply various compressive sensing algorithms
switch lower(params.csmethod) % [params] csmethod
    case 'gpsr'     % [2.1] gradient projection sparse reconstruction
        % [2.1.1] sparse representation
        rows = params.rows;
        cols = params.cols;
        switch lower(params.srbasis) % [params] srbasis
            case 'haar'    % Haar wavelet transform basis
                Fr = haarmtx(rows); % assist rows is power of 2
                if cols==rows
                    Fc = Fr;
                else
                    Fc = haarmtx(cols); % assist cols is power of 2
                end
                Psi = kron(Fc,Fr);  % sparsifing basis
            case 'dct'     % discrete cosine transform basis
                Dr = dctmtx(rows); % assist rows is power of 2
                if cols==rows
                    Dc = Dr;
                else
                    Dc = dctmtx(cols); % assist cols is power of 2
                end
                Psi = kron(Dc,Dr);  % sparsifing basis
            otherwise
                error('unsupported sparse representation basis %s.\n',params.srbasis);
        end
        A = double(sensmat)*Psi;
        
        % [2.1.2] ell_1 solver gpsr 
        % parameters configuration of gpsr
%         [M,N] = size(sensmat);
        sig0 = zeros(size(sensmat,2),1); % start point
        lambda = 1; % regularization parameter for ell_1-ell_2 optimization
%         sig0 = params.sig0; % start point
%         lambda = params.lambda; % regularization parameter for ell_1-ell_2 optimization
        
        tol = 1e-4; % convergence tolerance for gradient projection method
        opts.mu      = 0.1;   % [backtracking l.s.] (0,1/2)
        opts.beta    = 0.5;   % [backtracking l.s.] (0,1)
        opts.maxiter = 1000; % [iteration] maximum iteration
        % apply GPSR algorithm
        w = gpsr(A,meas,sig0,lambda,tol,opts);
        sig_out = Psi*w;
        % [end] GPSR
        %
        % 
    case 'tv'      % [2.2] total variation regularization (different from tval3)
        % [2.2.0] parameters configuration of tv
        rows = params.rows;
        cols = params.cols;
        % [2.2.1] generator tv operator (gradient operator)
        gradord = 1; % gradient order
        tvop = gentvop(rows,cols,gradord); 
        sensmat = double(sensmat);
        % [2.2.2] apply tv solver
        % opts.x0      = pinv(sensmat)*meas; % start point (pseudo-inverse result)
        opts.mu0     = 1;    % initial value of penalty factor mu
        opts.mubar   = 1e10; % maximum value of penalty factor mu
        opts.rho     = 1.05; % multiplication step of penalty factor mu
        opts.tol     = 1e-4;
        opts.miniter = 20;   % minimum iterations
        opts.maxiter = 500;  % maximum iterations
        % apply TV algorithm
        sig_out = tv(sensmat,meas,tvop,opts);
        % [end] TV
        % 
        % 
    case 'gap'     % [2.3] generalized alternating projection (GAP)
                   %       applying mixed Hadamard sensing method (MHSM)
        rows = params.rows;
        cols = params.cols;
        switch lower(params.srbasis) % [params] srbasis
            case 'haar'       % Haar wavelet transform basis
                Fr = haarmtx(rows); % assist rows is power of 2
                if cols==rows
                    Fc = Fr;
                else
                    Fc = haarmtx(cols); % assist cols is power of 2
                end
                Psi = kron(Fc,Fr);  % sparsifing basis
            case 'daubechies' % Daubechies wavelet transform basis (DB-8)
                level = 3;
                qmf = MakeONFilter('Daubechies',8); % Daubechies-8 wavelet
                sig_level_row = log2(rows);
                sig_level_col = log2(cols);
                Fr = get_waveletMatrix(qmf,sig_level_row,level,level);
                Fc = get_waveletMatrix(qmf,sig_level_col,level,level);
                Psi = kron(Fc,Fr);  % sparsifing basis
            case 'dct'        % discrete cosine transform basis
                Dr = dctmtx(rows); % assist rows is power of 2
                if cols==rows
                    Dc = Dr;
                else
                    Dc = dctmtx(cols); % assist cols is power of 2
                end
                Psi = kron(Dc,Dr);  % sparsifing basis
            otherwise
                error('unsupported sparse representation basis %s.\n',params.srbasis);
        end
        A = double(sensmat)*Psi;
        
        opts.tol     = 1e-3;
        opts.maxiter = 300;
        % apply GAP algorithm
        w = gap(A,meas,opts);
        sig_out = Psi*w;
        % [end] GAP
        % 
        % 
    case 'gap-tv'  % [2.4] generalized alternating projection (GAP) in 
                   % solving total variation (TV) minimization problem
        % [2.4.0] parameters configuration of tv
        rows  = params.rows;
        cols  = params.cols;
        % [2.4.1] generator tv operator (gradient operator)
        gradord = 1; % gradient order
        tvop = gentvop(rows,cols,gradord); 
        sensmat = double(sensmat);
        
        % [2.4.2] apply GAP-TV solver
        opts.lambda  = 16;        % TV regulizer
        opts.alpha   = 8;         % parameter for clipping algorithm, 
                                  %   alpha>=max(eig(tvop*tvop')) [=8]
        opts.tol     = 1e-5;      % convergence tolerent
        opts.acc     = false;     % accelerated version of GAP method
        % apply GAP_TV algorithm 
        sig_out = gap_tv(sensmat,meas,tvop,opts);
        % [end] GAP-TV
        % 
        % 
    case 'gap2d'   % [2.5] generalized alternating projection (GAP) for
                   %       two-dimensional (2D) compressive sensing
        % [2.5.1] expand the dimensions of rows and columns to be the power
        %         of two (2^k)
        rows = params.rows;
        cols = params.cols;
        rk   = ceil(log2(rows));
        ck   = ceil(log2(cols));
        rows_exp = 2^rk; % expanded number of rows (2^rk)
        cols_exp = 2^ck; % expanded number of columns (2^ck)
        sensmat = double(sensmat);
        [M,N] = size(sensmat);
        semsmat_exp = [sensmat zeros(M,rows_exp*cols_exp-N)];
        
        % [2.5.2] parameters configuration of GAP2D
        spbasis.space    = 'wavelet'; % transform for space, 'wavelet' or 'dct'
        spbasis.time     = 'dct';     % transform for spectrum, 'wavelet' or 'dct', dct is always used. I f we use wavelet, T need to be the power of 2. we can use no, means no transformation
        % spbasis.spectrum = 'no';  % Here we use no, means no transfromation in spectrum
        weighttype.space = 'tree';    % Here we can select:  'tree' or 'block'
        weighttype.time  = 'block';   % Here we can select:  'tree' or 'block', if we use tree, T need to be the power of 2
        weight_base.type = 'cos'; % here we can select 'exp' or 'cos'
        if strcmp(weight_base.type,'exp')
            weight_base.space = 1.5;   % This weight is the base of exponential decay. should be larger than 1 [1 2] is always used
            weight_base.time  = 1.5;
            weight_base.T     = 1.5;
        end
        % The block size for group
        block.row = 2;
        block.col = 2;
        % block.T   = T/2;
        block.T   = 1; % for 2D-CS
        % stop criterion
        stopc.iternum = 10^3;
        stopc.err     = 10^-5;
        acc           = 2; % GAP with acceleration
        % acc           = 0; % GAP without acceleration
        ydim          = rows_exp*cols_exp;
        m_star = ceil(ydim/(block.row*block.col*block.T));
        
        A_cs  = @(x) fm_cs(semsmat_exp,x);
        At_cs = @(y) fmt_cs(semsmat_exp,y);
        
        T = 1; % number of frames varying time
        m_star = m_star-1; % m_star should be decreased for 2D-CS
        % [2.5.3] apply GAP2D algorithm
        theta0 = gap2d(meas,A_cs,At_cs,rows_exp,cols_exp,T,block,...
                       spbasis,m_star,stopc,acc,weight_base,weighttype);
        
        theta_vec = theta0(:);
        sig_out = theta_vec(1:rows*cols); % recover the reconstructed signal
        % [end] GAP2D
        % 
        % 
    case 'gap3d'   % [2.6] generalized alternating projection (GAP) in 
                   %       three dimensions (3D)
        % 
        % to apear in CS3D
        % 
        % 
        % [end] GAP3D
        % 
        % 
    case 'ap'     % [2.7] alternating projection (AP) applied in 
                  %       single-pixel imaging
        % channel-wise recovery
        if isfield(params, 'x0')
            x0 = params.x0;
        end
        nch = size(meas, 2); % number of measurement channels
        for ich = 1:nch
            if isfield(params, 'x0')
                params.x0 = x0(:,ich);
            end
            sig_out(:,ich) = ap(sensmat, meas(:,ich), params); % note: params is transfered here
        end
        % [end] AP
        % 
        % 
    case 'tgi'    % [3.1] traditional ghost imaging (based on first-order correlation)
        [M,N] = size(sensmat);
        measmean = 0;
        pattmean = zeros(1,N);
        obj = pattmean;
        for i = 1:M
            measmean = (measmean*(i-1)+meas(i))/i;
            pattmean = (pattmean*(i-1)+double(sensmat(i,:)))/i;
            obj = (obj*(i-1)+(meas(i)-measmean)*(double(sensmat(i,:))-pattmean))/i;
        end
        sig_out = obj';
        % [end] TGI
        %
        %
    case 'zero-filling' % [4.1] Zero-filling for orthogonal bases
                        %       e.g., Walsh-Hadamard transform [-1/+1],
                        %       discrete cosine transform [-/+], and
                        %       discrete Fourier transform [-1/+1].
        sig_out = zero_filling(meas, params); % note: params is transfered here
        % [end] ZERO-FILLING
        % 
        % 
    case 'nlr-cs' % [5.1] Non-local mean regularization-based compressive 
                  %       sensing (NLR-CS)
        %%% TODO %%%
        % 
    case 'damp' % [5.2] Denoising-based approximate message passing (D-AMP)
                %       It works only with mixed Walsh-Hadamard transform, 
                %        Rademacher distributed patterns [-1/+1], and
                %        Gaussian distributed pattens [-/+] right now.
                %       It is also adapted to work with random binary
                %       pattners [0/1] ~ Bernoulli(1/2) after 
                %       mean-subtraction, which is sub-optimal and does not 
                %       yield comparable results as random [-/+] patterns 
                %       listed above.
        maxiter = params.maxiter;
        width   = params.cols;
        height  = params.rows;
        denoiser = params.denoiser;
        if min(sensmat(:))>=0 % random binary patterns
            sensmat = 2*sensmat-1;
            meas = 2*(meas-mean(meas));
        end
        alpha = sqrt(sum(abs(sensmat(:,1)).^2)); % normalization factor to the sensing matrix
        sensmat = sensmat/alpha;
        meas = meas/alpha;
        nch = size(meas, 2); % number of measurement channels
        % channel-wise reconstruction
        for ich = 1:nch
            sig_out(:,ich) = 1/255*DAMP(255*meas(:,ich), maxiter, height, width, denoiser, sensmat);
            % sig_out(:,ich) = DAMP(meas(:,ich), maxiter, height, width, denoiser, sensmat);
        end
        % [end] DAMP
        % 
        % 
    case 'pnp-ap' % [6.0] Plug-and-play algorithms based on 
                  %       alternating projection (PnP-AP)
        if isfield(params,'channel_wise') && params.channel_wise % [channel-wise]
            if isfield(params, 'x0')
                x0 = params.x0;
            end
            nch = size(meas, 2); % number of measurement channels
            for ich = 1:nch
                if isfield(params, 'x0')
                    params.x0 = x0(:,ich);
                end
                sig_out(:,ich) = pnp_ap(sensmat, meas(:,ich), params); % note: params is transfered here
            end
        else % [integrated]
            sig_out = pnp_ap(sensmat, meas, params); % note: params is transfered here
        end
        % [end] PnP-ADMM
        % 
        % 
    case 'pnp-admm' % [6.1] Plug-and-play algorithms based on 
                    %       alternating direction method of multipliers
                    %       (PnP-ADMM)
        if isfield(params,'channel_wise') && params.channel_wise % [channel-wise]
            if isfield(params, 'x0')
                x0 = params.x0;
            end
            nch = size(meas, 2); % number of measurement channels
            for ich = 1:nch
                if isfield(params, 'x0')
                    params.x0 = x0(:,ich);
                end
                sig_out(:,ich) = pnp_admm(sensmat, meas(:,ich), params); % note: params is transfered here
            end
        else % [integrated]
            sig_out = pnp_admm(sensmat, meas, params); % note: params is transfered here
        end
        % [end] PnP-ADMM
        % 
        % 
    case 'pnp-gap'  % [6.2] Plug-and-play algorithms based on generalized
                    %       alternating direction method (PnP-GAP)
        if isfield(params,'channel_wise') && params.channel_wise % [channel-wise]
            if isfield(params, 'x0')
                x0 = params.x0;
            end
            nch = size(meas, 2); % number of measurement channels
            for ich = 1:nch
                if isfield(params, 'x0')
                    params.x0 = x0(:,ich);
                end
                sig_out(:,ich) = pnp_gap(sensmat, meas(:,ich), params); % note: params is transfered here
            end
        else % [integrated]
            sig_out = pnp_gap(sensmat, meas, params); % note: params is transfered here
        end
        % [end] PnP-GAP
        % 
        % 
     case 'pnp-ist'  % [6.3] Plug-and-play algorithms based on iterative 
                     %       shrinkage-/thresholding (PnP-IST) algorithm
        if isfield(params,'channel_wise') && params.channel_wise % [channel-wise]
            if isfield(params, 'x0')
                x0 = params.x0;
            end
            nch = size(meas, 2); % number of measurement channels
            for ich = 1:nch
                if isfield(params, 'x0')
                    params.x0 = x0(:,ich);
                end
                sig_out(:,ich) = pnp_ist(sensmat, meas(:,ich), params); % note: params is transfered here
            end
        else % [integrated]
            sig_out = pnp_ist(sensmat, meas, params); % note: params is transfered here
        end
        % [end] PnP-IST
        % 
        % 
      case 'pnp-twist'  % [6.4] Plug-and-play algorithms based on two-step 
                        %       iterative shrinkage-/thresholding 
                        %       (PnP-TwIST) algorithm
        if isfield(params,'channel_wise') && params.channel_wise % [channel-wise]
            if isfield(params, 'x0')
                x0 = params.x0;
            end
            nch = size(meas, 2); % number of measurement channels
            for ich = 1:nch
                if isfield(params, 'x0')
                    params.x0 = x0(:,ich);
                end
                sig_out(:,ich) = pnp_twist(sensmat, meas(:,ich), params); % note: params is transfered here
            end
        else % [integrated]
            sig_out = pnp_twist(sensmat, meas, params); % note: params is transfered here
        end
        % [end] PnP-TwIST
        % 
        % 
      case 'pnp-fista'  % [6.5] Plug-and-play algorithms based on the fast 
                        %       iterative shrinkage-/thresholding algorithm
                        %       (PnP-FISTA) 
        if isfield(params,'channel_wise') && params.channel_wise % [channel-wise]
            if isfield(params, 'x0')
                x0 = params.x0;
            end
            nch = size(meas, 2); % number of measurement channels
            for ich = 1:nch
                if isfield(params, 'x0')
                    params.x0 = x0(:,ich);
                end
                sig_out(:,ich) = pnp_fista(sensmat, meas(:,ich), params); % note: params is transfered here
            end
        else % [integrated]
            sig_out = pnp_fista(sensmat, meas, params); % note: params is transfered here
        end
        % [end] PnP-FISTA
        % 
        % 
      case 'pnp-qcs' % [7.1] Plug-and-play algorithms based on 
                      %      Basis Pursuit DeQuantizer of mement p (BPDQ_p)
                      %      for quantized compressive sensing (QCS)
        if isfield(params,'channel_wise') && params.channel_wise % [channel-wise]
            if isfield(params, 'x0')
                x0 = params.x0;
            end
            nch = size(meas, 2); % number of measurement channels
            for ich = 1:nch
                if isfield(params, 'x0')
                    params.x0 = x0(:,ich);
                end
                sig_out(:,ich) = pnp_qcs(sensmat, meas(:,ich), params); % note: params is transfered here
            end
        else % [integrated]
            sig_out = pnp_qcs(sensmat, meas, params); % note: params is transfered here
        end
        % [end] PnP-QCS
        % 
        % 
    otherwise
        error('Unsupported compressive sensing algorithm %s.\n',params.csmethod);
end

