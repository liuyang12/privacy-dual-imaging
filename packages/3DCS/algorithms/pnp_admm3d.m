function [ x, psnrall ] = pnp_admm3d( A, y, opt )
%PNP_ADMM Plug-and-play (PnP) algorithm(s) based on alternating direction
%method of multipliers (ADMM) for compressive image reconstruction. 
% TODO: supposed to work with multi-channel measurements y's (and sensing
% matrices A's), e.g., color image with the same sensing matrix and video
% sequence with the same and different sensing matrices.
%   [x,psnrall]=PNP_ADMM(A,y,opt) returns the iterative solution for 
%   the ADMM solver and the PSNR in each iteration (if the ground truth 
%   opt.orig is available), where A is M-by-N sensing matrix, y is M-by-1 
%   measurement vector with (Gaussian) measurement noise, opt are options
%   of the solver attached with default settings. The expected range of x
%   is [0,1].
%   Note
%     For multi-channel measurements y's with size of M-by-F and the same
%     sensing matrix A, the size of A is still M-by-N (e.g.,color images),
%     whereas in multi-channel measurements with different sensing matrix
%     A, the size of A is M-by-N-by-F.
%   Model
%     Optimization model of the linear inverse problem
%        min  1/2*|A*x-y|_2^2+lambda*R(z),
%        s.t. x=z.
%     where |·|_2 is the l_2 norm for the data fidelity term
%     1/2*|A*x-y|_2^2, and R(·) is the regularization term, where
%     minimizing data fidelity and regularization terms individually 
%     would result in projection and denoising sub-problems, respectively.
%     And lambda is the balancing factor trading-off the regularization and
%     data fidelity. 
%     PnP-ADMM solution
%       x^{k+1} = (rho*I+A'*A)^{-1}*(A'*y+rho*(z^k-u^k))
%       z^{k+1} = D_{sigma_k}(x^{k+1}+u^k), sigma_k=lambda/rho
%       u^{k+1} = u^k + (x^{k+1}-z^{k+1})
%   Reference(s)
%     [2] S. Boyd, N. Parikh, E. Chu, B. Peleato, and J. Eckstein, 
%         Distributed Optimization and Statistical Learning via the 
%         Alternating Direction Method of Multipliers, Foundations and 
%         Trends® in Machine Learning, vol. 3, no. 1, pp. 1-122, 2011.
%   Copyright(C) 2019-2020 by <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>.
%   Last modified May 23, 2020.

if nargin<2
    opt = [];
end

% [0] default parameter configuration, to be specified
denoiser = 'tv';    % denoiser as PnP regularizer
% x0       = A'*(y); % start point (initialization of iteration)
lambda   = 1;       % balancing regularization factor
rho      = 1e-3;    % Lagrangian multiplier (balancing noise)
maxiter  = 100;     % maximum number of iteration
sigma    = 10/255;  % noise deviation 
nosestim = true;    % enable noise estimation (if possible)
tvweight = 0.07;    % weight for TV denoising
tviter   = 5;       % number of iteration for TV denoising
flag_iqa = true;    % flag of showing image quality assessments
ffdnetvnorm_init = true;  % use normalized image as input for initialization
                          %  (with the first 10 iterations)
bm3d_profile = 'np';      % quailty & complexity trade-off profile selection 
                          %  ('np' for normal profile, 'lc' for low-complexity)
orthogonal_basis = false; % sensing matrix A orthogonal or not 

if isfield(opt,'denoiser'), denoiser = opt.denoiser; end
if isfield(opt,'x0'),             x0 = opt.x0;       end
if isfield(opt,'lambda'),     lambda = opt.lambda;   end
if isfield(opt,'rho'),           rho = opt.rho;      end
if isfield(opt,'maxiter'),   maxiter = opt.maxiter;  end
if isfield(opt,'tvweight'), tvweight = opt.tvweight; end
if isfield(opt,'tviter'),     tviter = opt.tviter;   end
if isfield(opt,'nosestim'), nosestim = opt.nosestim; end
if isfield(opt,'sigma'),       sigma = opt.sigma;    end
if isfield(opt,'flag_iqa'), flag_iqa = opt.flag_iqa; end
if isfield(opt,'ffdnetvnorm_init'), ffdnetvnorm_init = opt.ffdnetvnorm_init; end
if isfield(opt,'bm3d_profile'), bm3d_profile = opt.bm3d_profile; end
if isfield(opt,'orthogonal_basis'), orthogonal_basis = opt.orthogonal_basis; end

if isfield(opt, 'nframe'), F = opt.nframe; else F = size(y,2); end
if ndims(A) <= 2, SAME_SENSING_MATRIX=true; else SAME_SENSING_MATRIX=false; end
if ~exist('x0','var') || isempty(x0) % start point
    for i = 1:F
        if SAME_SENSING_MATRIX % same sensing matrix for all measurements
            x0(:,i) = A'*y(:,i); 
        else
            x0(:,i) = A(:,:,i)'*y(:,i); 
        end
    end
end

if  isfield(opt,'ffdnetvnorm') && ffdnetvnorm_init % 
    sigma = [50/255 sigma];
    maxiter = [10 maxiter];
    ffdnetvnorm = opt.ffdnetvnorm;
end

% size of the image or the image sequence (used when resizing)
N = size(A,2);
n = sqrt(N);
if ~isfield(opt, 'rows') && ~isfield(opt, 'cols')
    imsize = [n n];
elseif isfield(opt, 'rows')
    imsize = [opt.rows N/opt.rows];
else
    imsize = [N/opt.cols opt.cols];
end
nch = size(y,2); % number of channels in the measurement
if nch > 1
    imsize = [imsize nch];
end

% pre-calculation
if orthogonal_basis % [orthogonal basis] A*A' diagonal
    if SAME_SENSING_MATRIX
        d = diag(A*A');
    else
        for i = 1:F
            d(:,i) = diag(A(:,:,i)*A(:,:,i)');
        end
    end
else % [general basis]
    I = eye(N); % identity matrix (N-by-N)
    if SAME_SENSING_MATRIX
        invmat = (rho*I+A'*A)\I; % (rho*I+A'*A)^{-1}
        x_rho = invmat*A'*y; % (rho*I+A'*A)^{-1}*A'*b
    else
        for i = 1:F
            invmat(:,:,i) = (rho*I+A(:,:,i)'*A(:,:,i))\I; % (rho*I+A'*A)^{-1}
            x_rho(:,i) = invmat(:,:,i)*A(:,:,i)'*y(:,i); % (rho*I+A'*A)^{-1}*A'*b
        end
    end
end

u  = zeros(size(x0),'like',x0); % residual

% [1] start iteration
z = x0; % auxiliary variable [initialization]
psnrall = []; % return empty with no ground truth
k = 1; % current number of iteration
for isig = 1:length(maxiter) % extension for a series of noise levels
    nsigma = sigma(isig); 
    opt.sigma = nsigma;
    for iter = 1:maxiter(isig)
        % [1.1] Euclidean projection
        if orthogonal_basis % [orthogonal basis] A*A' diagonal
            if SAME_SENSING_MATRIX
                x = (z-u) + A'*( (y-A*(z-u))./(rho+d) ); % element-wise
            else
                for i = 1:F
                    x(:,i) = (z(:,i)-u(:,i)) + A(:,:,i)'* ( (y(:,i)-A(:,:,i)*(z(:,i)-u(:,i)))./(rho+d(:,i)) );
                end
            end
        else % [general basis]
            if SAME_SENSING_MATRIX
                x = x_rho + rho*invmat*(z-u); % matrix product
            else
                for i = 1:F
                    x(:,i) = x_rho(:,i) + rho*invmat(:,:,i)*(z(:,i)-u(:,i));
                end
            end
        end
        xu_mat = reshape(x+u, imsize); % reshape to image format
        if ~isreal(xu_mat)
            xu_mat = real(xu_mat);
            % xu_mat = abs(xu_mat);
        end
        % [1.2] Denoising to match the video prior
        switch lower(denoiser)
            case 'tv' % TV denoising
                z_mat = TV_denoising(xu_mat,tvweight,tviter);
            case 'bm3d' % BM3D denoising
                [~,z_mat] = BM3D(1,xu_mat*255,nsigma*255,bm3d_profile,0); % nsigma
            case 'cbm3d' % BM3D denoising
                [~,z_mat] = CBM3D(1,xu_mat*255,nsigma*255,bm3d_profile,0); % nsigma
            % case 'wnnm' % WNNM image denoising (MATLAB-style matrix version)
            %     z_mat = wnnm_imdenoise(xu_mat,[],opt); % opt.sigma
            % case 'ffdnet' % FFDNet image denoising
            %     if ffdnetnorm_init
            %         if isig==1
            %             opt.ffdnetnorm = true;
            %         else
            %             opt.ffdnetnorm = ffdnetnorm;
            %         end
            %     end
            %     z_mat = ffdnet_imdenoise(xu_mat,[],opt); % opt.sigma
                
            case 'vbm3d' % VBM3D denoising
                [~,z_mat] = VBM3D(xu_mat,nsigma,0,0);
            case 'vbm4d' % VBM4D denoising
                if nosestim % noise estimation enabled
                    z_mat = vbm4d(xu_mat,-1,'lc',1,1,1,0); % -1 to enable noise estimation
                else % noise estimation disabled
                    z_mat = vbm4d(xu_mat,nsigma,'lc',1,1,1,0); % sigma as estimated noise level
                end
            case 'bm4d' % BM4D denoising
                if nosestim % noise estimation enabled
                    z_mat = bm4d(xu_mat,'Gauss',0,'lc',1,0); % 0 to enable noise estimation
                else % noise estimation disabled
                    z_mat = bm4d(xu_mat,'Gauss',nsigma,'lc',1,0); % sigma as estimated noise level
                end
            case 'wnnm' % WNNM video denoising 
                z_mat = wnnm_vdenoise(xu_mat,[],opt); % opt.sigma as estimated noise level
            case 'ffdnet' % FFDNet video denoising (frame-wise)
                if ffdnetvnorm_init
                    if isig==1
                        opt.ffdnetvnorm = true;
                    else
                        opt.ffdnetvnorm = ffdnetvnorm;
                    end
                end
                z_mat = ffdnet_vdenoise(xu_mat,[],opt); % opt.sigma
            otherwise
                error('Unsupported denoiser %s!',denoiser);
        end
        z = reshape(z_mat, size(x0)); % reshape back to vector format
        u = u + (x-z); % update residual
        % [1.3] save and show intermediate results of psnr and ssim
        if flag_iqa && isfield(opt,'orig') && (~isempty(opt.orig))
            psnrall(k) = psnr(double(x),double(opt.orig)); % record all psnr
            % ssimall(k) = ssim(double(v),opt.orig); % record all ssim
            if (mod(k,5)==0) 
                fprintf('  PnP-ADMM-%s iteration % 4d, sigma % 3d, PSNR %2.2f dB.\n',...
                    upper(opt.denoiser),k,nsigma*255,psnrall(k));
            end
        end
        k = k+1;
    end % ADMM loop [maxiter]
end % sigma loop [length(maxiter)]

end            
            
            