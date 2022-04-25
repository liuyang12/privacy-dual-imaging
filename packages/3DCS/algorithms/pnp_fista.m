function [ x, psnrall ] = pnp_fista( A, y, opt )
%PNP_FISTA Plug-and-play (PnP) algorithm(s) based on the fast iteraive 
%shrinkage-thresholding algorithm (FISTA) for compressive image 
%reconstruction. 
% TODO: supposed to work with multi-channel measurements y's (and sensing
% matrices A's), e.g., color image with the same sensing matrix and video
% sequence with the same and different sensing matrices.
%   [x,psnrall]=PNP_TWIST(A,y,opt) returns the iterative solution for 
%   the FISTA solver and the PSNR in each iteration (if the ground truth 
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
%        min  1/2|A*x-y|_2^2 + lambda*R(x)
%     where A*x=y is the linear manifold for the data fidelity term, and 
%     R(·) is the regularization term, where minimizing data fidelity and 
%     regularization terms individually would result in projection and 
%     denoising sub-problems, respectively.
%     PnP-FISTA solution
%       x^{k+1} = theta^k + 1/L*A'*(y-A*theta^{k})
%       z^{k+1} = D_{sigma_k}(x^{k+1})
%       theta^{k+1} = z^{k+1} + (t^k-1)/t^{k+1}*(z^{k+1}-z^k)
%       where 
%         t^{k+1}=(1+sqrt(1+4*{t^k}^2))/2 and t^1=1.
%   Reference(s)
%     [1] A. Beck and M. Teboulle, A Fast Iterative Shrinkage-Thresholding 
%         Algorithm for Linear Inverse Problems, SIAM Journal on Imaging 
%         Sciences, vol. 2, no. 1, pp. 183-202, 2009.

%   Copyright(C) 2019-2020 by <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>.
%   Last modified July 3, 2020.

if nargin<2
    opt = [];
end

% [0] default parameter configuration, to be specified
denoiser = 'tv';    % denoiser as PnP regularizer
% x0       = A'*(y); % start point (initialization of iteration)
lambda   = 1;       % balancing regularization factor
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
if isfield(opt,'maxiter'),   maxiter = opt.maxiter;  end
if isfield(opt,'tvweight'), tvweight = opt.tvweight; end
if isfield(opt,'tviter'),     tviter = opt.tviter;   end
if isfield(opt,'nosestim'), nosestim = opt.nosestim; end
if isfield(opt,'sigma'),       sigma = opt.sigma;    end
if isfield(opt,'flag_iqa'), flag_iqa = opt.flag_iqa; end
if isfield(opt,'ffdnetvnorm_init'), ffdnetvnorm_init = opt.ffdnetvnorm_init; end
if isfield(opt,'bm3d_profile'), bm3d_profile = opt.bm3d_profile; end
if isfield(opt,'orthogonal_basis'), orthogonal_basis = opt.orthogonal_basis; end

if ~exist('x0','var') || isempty(x0) % start point
    x0 = A'*y; 
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
% Lipschitz constant of gradient of f (1/2*|A*x-y|_2^2) 
if orthogonal_basis % [orthogonal basis] A*A' diagonal
    L = max(diag(A*A')); % Lipschitz constant of gradient of f (1/2*|A*x-y|_2^2) 
elseif size(A,1)<=size(A,2) % [general basis] sub-sampled or complete-sampled
    L = max(eig(A*A')); % Lipschitz constant of gradient of f (1/2*|A*x-y|_2^2) 
else
    L = max(eig(A'*A)); % Lipschitz constant of gradient of f (1/2*|A*x-y|_2^2) 
end

% [1] start iteration
theta = x0; % auxiliary variable [initialization]
z = x0; % previous value of theta
t = 1;
psnrall = []; % return empty with no ground truth
k = 1; % current number of iteration
for isig = 1:length(maxiter) % extension for a series of noise levels
    nsigma = sigma(isig); 
    opt.sigma = nsigma;
    for iter = 1:maxiter(isig)
        % [1.1] Euclidean projection
        yb = A*theta;
        x = theta + 1/L*A'*(y-yb); % IST with fixed step size
        
        z_prev = z;
        t_prev = t;
        
        x_mat = reshape(x, imsize); % reshape to image format
        % [1.2] Denoising to match the video prior
        switch lower(denoiser)
            case 'tv' % TV denoising
                z_mat = TV_denoising(x_mat,tvweight,tviter);
            case 'bm3d' % BM3D denoising
                [~,z_mat] = BM3D(1,255*x_mat,nsigma*255,bm3d_profile,0); % nsigma
            case 'cbm3d' % BM3D denoising
                [~,z_mat] = CBM3D(1,255*x_mat,nsigma*255,bm3d_profile,0); % nsigma
            case 'wnnm' % WNNM video denoising (MATLAB-style matrix version)
                z_mat = wnnm_imdenoise(x_mat,[],opt); % opt.sigma
            case 'ffdnet' % FFDNet video denoising (frame-wise)
                if ffdnetvnorm_init
                    if isig==1
                        opt.ffdnetvnorm = true;
                    else
                        opt.ffdnetvnorm = ffdnetvnorm;
                    end
                end
                z_mat = ffdnet_imdenoise(x_mat,[],opt); % opt.sigma
            otherwise
                error('Unsupported denoiser %s!',denoiser);
        end
        z = reshape(z_mat, size(x0)); % reshape back to vector format
        t = (1+sqrt(1+4*t^2))/2;
        theta = z + (t_prev-1)/t * (z-z_prev); % FISTA update
        % [1.3] save and show intermediate results of psnr and ssim
        if flag_iqa && isfield(opt,'orig') && (~isempty(opt.orig))
            psnrall(k) = psnr(double(x),double(opt.orig)); % record all psnr
            % ssimall(k) = ssim(double(v),opt.orig); % record all ssim
            if (mod(k,5)==0) 
                fprintf('  PnP-FISTA-%s iteration % 4d, sigma % 3d, PSNR %2.2f dB.\n',...
                    upper(opt.denoiser),k,nsigma*255,psnrall(k));
            end
        end
        k = k+1;
    end % FISTA loop [maxiter]
end % sigma loop [length(maxiter)]

end            
            
            