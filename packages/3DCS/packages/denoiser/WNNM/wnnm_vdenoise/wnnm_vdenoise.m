function [denoisedv,psnrall,ssimall] = wnnm_vdenoise(noisyv,orgv,para)
%WNNM_VDENOISE Weighted nuclear norm minimization (WNNM)-based video
%denoising.
%   denoisedv=WNNM_VDENOISE(noisyv,para) returns the denoised image
%   denoisedv using WNNM, where noisyv is the input noisy video and para
%   is the parameters for the WNNM video denoiser.
%   # Model
%     The basic idea of WNNM for video denoising is that the patch group
%     with nonlocal self-similarity can be expressed as a low-rank matrix.
%     Low-rank matrix approximation problem has already been solved by
%     low-rank matrix factorization (LRMF) and nuclear norm minimization 
%     (NNM). WNNM approach to this problem focuses on varied weights of
%     each singular value, where greater weight indicates less importance.
%     The CVPR'14 paper gives an convergence proof to the fixed point of 
%     the WNNM problem with non-descending weights, which is common in 
%     image denoising applications where small singular values of the patch
%     group are supposed to be surpressed.
%   # Parameters
%     parameters              [default] value
%       para.patchsize
%   # References
%   [1]  S. Gu, L. Zhang, W. Zuo, and X. Feng, "Weighted Nuclear Norm 
%          Minimization with Application to Image Denoising," in 2014 IEEE 
%          Conference on Computer Vision and Pattern Recognition (CVPR), 
%          2014, pp. 2862-2869.
%   [2]  S. Gu, Q. Xie, D. Meng, W. Zuo, X. Feng, and L. Zhang, "Weighted 
%          Nuclear Norm Minimization and Its Applications to Low Level 
%          Vision," International Journal of Computer Vision, vol. 121, 
%          no. 2, pp. 183-208, 2017.
%   See also .

%   Code copyright(C) 2017-2018 by <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>. This code is a
%   derivative of WNNM code <a href="matlab: 
%   web('http://www4.comp.polyu.edu.hk/~cslzhang/code/WNNM_code.zip')"> 
%   provided by S. Gu. Note that the original WNNM code is for image
%   denoising only, and this code is a non-trivial attempt for video
%   denoising.
%   Last modified Dec 12, 2017.

if isfield(para,'enparfor') && para.enparfor % enable parfor
    para.enparfor = true;
    [denoisedv,psnrall,ssimall] = wnnm_vdenoise_par(noisyv,orgv,para);
    return;
end

para = vdefparaconf(para); % default parameter configuration

[nrow,ncol,nframe] = size(noisyv); % grayscale video (color todo)
% [1] pre-calculation of the indexes of the neighbors within the search
%     window
[neighborindarr,neighbornumarr,selfindarr] = neighborind([nrow ncol],para);

% [2] WNNM denoisng for several iterations
estv = noisyv;
psnrall = [];
ssimall = [];
for iter = 1:para.iternum
    % correction between adjacent iterations
    if ~para.adaptboost % 
        estv = estv + para.delta*(noisyv-estv); 
    else % adaptive boosting for WNNM-based denoising
        Eestv = mean(abs(estv(:)).^2);
        Enoise  = para.abeta*abs(para.nsigma^2-var(noisyv(:)-estv(:)));
        rho = sqrt(Eestv)/(sqrt(Eestv)+sqrt(max(Eestv-Enoise,0)));
        fprintf('    Iteration % 2d, rho = %.3f.\n',iter,rho);
        estv = estv + (1-rho)*(noisyv-estv); 
    end
    v = zeros(nrow,ncol,nframe);
    f = zeros(nrow,ncol,nframe);
    for iframe = 1:nframe
        lframe = max(1,min(max(1,iframe-para.sframesize),...
            nframe-2*para.sframesize));
        rframe = min(nframe,max(min(nframe,iframe+para.sframesize),...
            2*para.sframesize+1));
        curseq = lframe:rframe;
        % current video excerpt to be denoised
        curnosv = noisyv(:,:,curseq); 
        curestv = estv(:,:,curseq);
        % use all frames for denoising
        cframe = iframe-lframe+1;
        
        % [2.1] splite the whole image into overlapped patches
        [rawpatchmat,nsigmamat] = v2patch(curestv,curnosv,para); 
        % inner loop to reuse the block-matching results
        
            if iter==1 % initial noise level of each patch
                nsigmamat = para.nsigma*ones(size(nsigmamat));
            end
            
        % [2.2] calculate the with non-local similarity for each 
        %       key patch
        [nonlocalarr,vindarr] = vblockmatch(cframe,rawpatchmat,...
            neighborindarr,neighbornumarr,selfindarr,para);
        % [2.3] patch estimation by means of WNNM
        [estpatchmat,frqpatchmat] = vpatchestimate(cframe,nonlocalarr,...
            vindarr,rawpatchmat,nsigmamat,selfindarr,para);
        % [2.4] aggregate overlapped patches to the whole image
        [curv,curf] = patch2v(estpatchmat,frqpatchmat,size(noisyv),...
            para.patchsize);
        v(:,:,curseq) = v(:,:,curseq)+curv;
        f(:,:,curseq) = f(:,:,curseq)+curf;
    end
    
    estv = v./(f+eps);
    
    if mod(iter-1,para.innerloop)==0 
        % % [2.2] calculate the patches with non-local similarity for 
        % %       each key patch
        % [nonlocalarr,vindarr] = vblockmatch(cframe,rawpatchmat,...
        %     neighborindarr,neighbornumarr,selfindarr,para);
        % less non-local patches with lower noise level
        para.patchnum = para.patchnum-10; 
    end
    if ~isempty(orgv)
        % error statistics and reconstruction quality accessment
        psnrall(iter) = vpsnr(estv,orgv,255);
        ssimall(iter) = vssim(estv,orgv);
        fprintf('  WNNM iteration % 2d, mean PSNR %2.2f dB, mean SSIM %.4f.\n',...
            iter,psnrall(iter),ssimall(iter));
    end
end
denoisedv = estv;

end

