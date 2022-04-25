% get the ghost image
% input:
%   - P: H x W x N
%   - S: N x 1, sample signal
%   - R: N x 1, reference signal
% output:
%   - ghost image
% v3: 1. add a new CS method CS_BM3D. 2. seprate TVAL3 from old CS method.
%     3. add a new GI method GSGI(Gerchberg-Saxton-like ghost imaging).

function [O, O_tval3] = cs_image(P,S,R,params)


[im_h,im_w,~] = size(P);
patterns_num = length(S);

O = zeros(im_h,im_w);
O_tval3 = zeros(im_h,im_w);

%% GI defination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R_sum = zeros(im_h,im_w);
% B_sum = 0;
% BI_sum = zeros(im_h,im_w);
%
% for i = 1:patterns_num
%
%     im = P(:,:,i);
%
%     B_sum = B_sum + S(i);
%     R_sum = R_sum + im;
%     BI_sum = BI_sum + im*S(i);
% end
%
% O = 1/patterns_num*BI_sum - 1/patterns_num*B_sum*1/patterns_num*R_sum;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch params.GI_method
    %% TGI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'TGI'
        S_mean = 0;
        P_mean = zeros(im_h,im_w);
        
        for i = 1:patterns_num
            
            S_mean = (S_mean*(i-1) + S(i))/i;
            P_mean = (P_mean*(i-1) + P(:,:,i))/i;
            
            O = (O*(i-1) + (S(i) - S_mean)*(P(:,:,i) - P_mean))/i;
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% DGI
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'DGI'
        S_mean = 0;
        R_mean = 0;
        P_mean = zeros(im_h,im_w);
        
        for i = 1:patterns_num
            
            S_mean = (S_mean*(i-1) + S(i))/i;
            P_mean = (P_mean*(i-1) + P(:,:,i))/i;
            R_mean = (R_mean*(i-1) + R(i))/i;
            
            O = (O*(i-1) + (S(i) - S_mean/R_mean*R(i))*(P(:,:,i) - P_mean))/i;
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% NGI
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'NGI'
        S_mean = 0;
        P_mean = zeros(im_h,im_w);
        R_mean = 0;
        
        for i = 1:patterns_num
            
            S_mean = (S_mean*(i-1) + S(i))/i;
            P_mean = (P_mean*(i-1) + P(:,:,i))/i;
            R_mean = (R_mean*(i-1) + R(i))/i;
            
            O = (O*(i-1) + (S(i)/R(i) - S_mean/R_mean)*(P(:,:,i) - P_mean))/i;
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'GSGI' 
        %%% Gerchberg-Saxton-like ghost imaging %%%
        % parameters for GSGI
        if ~isfield(params, 'err_threshold')
            params.err_threshold = 1e-10;   % mean-square error threshold of the reconstructed Obj.
        end
        if ~isfield(params, 'max_iter');
            params.max_iter = 500;  % maxmium iteration times
        end
        S_mean = 0;
        P_mean = zeros(im_h,im_w);
        
        for i = 1:patterns_num
            
            S_mean = (S_mean*(i-1) + S(i))/i;
            P_mean = (P_mean*(i-1) + P(:,:,i))/i;
            
            O = (O*(i-1) + (S(i) - S_mean)*(P(:,:,i) - P_mean))/i; % TGI result as initial value
            
        end
%         O = O_first;
%         O = O*S_mean;
%         O = rand(im_h, im_w);
        % normalization
%         O = (O - min(min(O)))/(max(max(O)) - min(min(O)));
        O_init = O;
        end_flag = 0;   % flag of iteration converges
        k = 1;
        for iter = 1:params.max_iter
            if end_flag
                break;
            end
            for m = 1:patterns_num
                Pm = P(:,:,m);
                Ptm = Pm.*O;
                Ftm = fft2(Ptm);
                % compensation of Ftm
% %                 Ftm_comp = Ftm - Ftm(1,1) + S(m);
%                 Ftm_comp = Ftm;
%                 Ftm_comp(1,1) = S(m);
%                 Ptm_upd = ifft2(Ftm_comp);
                Ptm_upd = Ptm + (S(m)-sum(sum(Ptm)))/(im_h*im_w);
                O_upd = O + Pm.*(Ptm_upd - O.*Pm)/max(max(max(abs(Pm))))^2;
                errall(k,1) = sum(sum((O_upd - O).^2))/(im_h*im_w);
                k = k+1;
                O = O_upd;
            end
            % mean-square error of updated and original Obj. 
%             errall(iter,1) = sum(sum((O_upd - O).^2))/(im_h*im_w);
            if sum(sum((O_upd - O).^2))/(im_h*im_w) <= params.err_threshold
%                 end_flag = 1;
%                     break;
            end
            iter
        end
        figure; plot(errall);
        O = O;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'APGI' % alternating projection ghost imaging
        O = fun_SPI_R_AP(P,S);
        % [case end] APGI
    case 'TVAL3'
        opts=[];opts.nonneg=true;opts.TVnorm = 1;
        opts.scale_A=false; opts.scale_b=false;
        
        if isfield(params, 'TVAL3_mu')
            opts.mu = params.TVAL3_mu;
        end
        if isfield(params, 'TVAL3_beta')
            opts.beta = params.TVAL3_beta;
        end
        
        row = size(P,1);
        col = size(P,2);
        P_vec = reshape(P,[row*col, size(P,3)])';
        [O,~]=TVAL3( P_vec,S,row,col,opts);
        %%%%%%%%%%%%% tune the parameter of mu and beta of TVAL3 %%%%%%%%
%          for k = 6:16
%             opts.mu = 2^k;
%             for kk = 1:7
%                 opts.beta = 2^kk;
%                 [img_tval3,~]=TVAL3( P,z,row,col,opts);
%                 O_tval3 = img_tval3;
%                 x_tval3=img_tval3(:);
%                 
%                 O_tval3 = (O_tval3 - min(min(O_tval3)))/((max(max(O_tval3))) - min(min(O_tval3)));
%                 O_tval3 = imresize(O_tval3,10,'nearest');
%                 
%                 mkdir(params.out_path);
%                 imwrite(im2uint8(O_tval3),sprintf('%s/O_tval3_mu%08d_beta%05d.png',params.out_path,opts.mu,opts.beta));
%             end
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'PINV'
        row = size(P,1);
        col = size(P,2);
        P_vec = reshape(P,[row*col, size(P,3)])';
        
        O = pinv(double(P_vec))*S;
        O=reshape(O,[row,col]);
    case 'CS'
        % parameters
        %         mu_0 = 1;
        %         mu_bar = 1e7;
        %         rho = 1.05;
        %         max_iter = 500;
        %         tol = 1e-3;
        %         plot_flag = 0;
        
        mu_0 = 1;
        mu_bar = 1e7;
        rho = 1.0;
        max_iter = 500;
        tol = 1e-9;
        plot_flag = 1;
        
        %%% load in pattern -------------------------------------------------------
        row=size(P,1);
        col=size(P,2);
        pattern_number=size(P,3);
        pattern=reshape(P,[row*col,pattern_number]);
        pattern=double(pattern');
        
        % the operator for calculating gradient
        gradient_order = 1;
        psi_tv = generate_gradient_operator(row,col,gradient_order);
        
        % measurement number
        m=pattern_number;
        P=pattern(1:m,:);
        
        % load in real data captured data
        z = S;% here input your captured data
        
        % get the initial solution with the algorithm: TVAL3
        opts=[];opts.nonneg=true;opts.TVnorm = 1;
       
     
        
        disp('calculating pinv_mat');
        pinv_mat_tv=pinv(psi_tv'*psi_tv+P'*P);
        
        disp('compressive sensing');
        [x_tv]=ccgi(psi_tv,P,z,[],mu_0,rho,tol,20,...
            plot_flag,row,col,pinv_mat_tv);
        
        img_tv=reshape(x_tv,[row,col]);
        O = img_tv;
    case 'CS_v2'
        % parameters
        mu_0 = 1; 
        mu_bar = 1e7;
        rho = 1.05;
        max_iter = 500;
        tol = 1e-3;
        plot_flag = 0;
        
        %%% load in pattern -------------------------------------------------------
        row=size(P,1);
        col=size(P,2);
        pattern_number=size(P,3);
        pattern=reshape(P,[row*col,pattern_number]);
        pattern=double(pattern');
        
        % the operator for calculating gradient
        gradient_order = 1;
        psi_tv = generate_gradient_operator(row,col,gradient_order);
        
        % measurement number
        m=pattern_number;
        P=pattern(1:m,:);
        
        % load in real data captured data
        z = S;% here input your captured data
        
        % get the initial solution with the algorithm: TVAL3
        opts=[];opts.nonneg=true;opts.TVnorm = 1;
        [img_tval3,~]=TVAL3( P,z,row,col,opts);
        x_tval3=img_tval3(:);
        
        disp('calculating pinv_mat');
        pinv_mat_tv=pinv(psi_tv'*psi_tv+P'*P);
        
        disp('compressive sensing');
        [x_tv]=ccgi(psi_tv,P,z,x_tval3,mu_0,rho,tol,20,...
            plot_flag,row,col,pinv_mat_tv);
        img_tv=reshape(x_tv,[row,col]);
        O = img_tv;
        
    case 'CS_BM3D'
        
        yN = size(P,1);
        xN = size(P,2);
        K = size(P,3);
        
        
        B = double(reshape(P,[xN*yN,K]));
        zz = S;
        
        % all parameters
        IterNumber=100;
        
        d_factor=3; mu=2;
        threshold_ampl=.05; sigma=0.; gamma_00=.2; % PSNR=45.82 dB
        
        Nstep = 1;  N11=4;  N22=8;
        
        threshType='h';  %% 'h' for hard-thresholding; 's' for soft-thresholding % hard-thresholding is recommended
        filtering=1;       %% filtering=1 switches on BM3D filtering; filtering=0 switches off BM3D filtering
        modulation_type=1;
        
        KAPPA=1000000;
        
        Gauss_var=0;            %% '1' for the varying weights and '0' for the invariant weights of approximative Gaussian distribution
        
        IdealWeights=1;        %% Ideal varying weights for Gaussian distribution with the varying variances; THIS PARAMETER WORKS only for Gauss_var=1;
        IdealWeights=0;         %% Adaptive weights for var GAUSS; THIS PARAMETER WORKS only for Gauss_var=1
        
        ss00=floor(IterNumber/4);  %% update FRAMES and go to varying Gaussian estimation
        
        
        ssigma_2=(mean(zz(:))/KAPPA);  % estimate of variance
        gamma_0=gamma_00/ssigma_2;  % parameter o the algorithm
        
        FHI=pinv(B'*B+eye(K)/gamma_0);
        delta_x_zz=FHI*(B'*B)*zz;
        
        mean_B=mean(B');   %% for sigma^2 update
        
        %% Approximative Gaussian ALGORITHM, CSGI_inv and CSGI_var
        
        v0=ones(yN,xN)/2;          %% Intialization
        
        
        for ss= 1:IterNumber                           %% Algorithm's iterations
            
            fprintf('%d. ',ss);
            if mod(ss,10)  == 0
                fprintf('\n');
            end
            
            if (Gauss_var==1)&(ss>=ss00)           %% START CSGI-varying (variang variance)
                if (ss==ss00)&(IdealWeights==1)
%                     'Ideal Weights'
                    D=diag(1./z_noiseless*KAPPA)*ssigma_2;
                    FHI=pinv(B'*B*D+eye(K)/gamma_0);
                    
                    delta_x_zz=FHI*(B'*B)*D*zz;
                end
                %% CSGI_var
                if (ss>=ss00)&(IdealWeights==0)
                    
                    for ss1=1:K
                        z_noiseless(ss1)=sum(B(:,ss1).*abs(u0_est(:)));
                    end
                    D=diag(1./z_noiseless*KAPPA)*ssigma_2;
                    FHI=pinv(B'*B*D+eye(K)/gamma_0);
                    
                    delta_x_zz=FHI*(B'*B)*D*zz;
                end
                %% end of CSGI_var
                
                x=delta_x_zz+FHI*B'*v0(:)/gamma_0;
                
                u0_est=v0(:)+B*D*(zz-x)*gamma_0;
                
                
            else
                %% CSGI_inv s
                x=delta_x_zz+FHI*B'*v0(:)/gamma_0;
                
                u0_est=v0(:)+B*(zz-x)*gamma_0; %0*PSI_PSI*rand(yN*yN,1)/(ss);
                
                %% variance update ( can be dropped)
                ssigma_22(ss)=mean(B')*abs(u0_est)/KAPPA;
                gamma_0=gamma_00/ssigma_2;
                
            end
            
            u0_est=reshape(u0_est,[yN,xN]);
            u0_est=abs(u0_est).*double((real(u0_est)>=0));
            
            %% RESULTS of FIRST STEP
%             ee_abs=u0_est-abs(u0);            
%             rmse_abs(ss)=sqrt(mean(ee_abs(:).^2));
            
            %% BM3D-frame filtering
            %keyboard
            %% FILTERING STARS from ITERATION ss=1
            if (ss>=1)&(ss<ss00)
                
                kk=10;  u0_est_filt=padarray((u0_est), [kk kk],'replicate', 'both'); %% PADDING FOR FILTERING, can be k=0, NO PADDING;
                
                [v0tmp1] = BM3D_filter(u0_est_filt,threshType, threshold_ampl, N11, N22, Nstep,filtering,ss,1);
                
                v0=v0tmp1(kk+1:end-kk,kk+1:end-kk);
                
            else
                %% Frame update for ss=ss00
                kk=10;  u0_est_filt=padarray((u0_est), [kk kk],'replicate', 'both'); %% PADDING FOR FILTERING, can be k=0, NO PADDING;
                
                [v0tmp1] = BM3D_filter(u0_est_filt,threshType, threshold_ampl, N11, N22, Nstep,filtering,ss,ss00);
                
                
                v0=v0tmp1(kk+1:end-kk,kk+1:end-kk);
                
                
                % u0_est=v0;
            end
            
            %% ITERATION VISUALIZATION
            O_show = imresize(abs(u0_est),10,'nearest');
            O_show = (O_show - min(min(O_show)))/(max(max(O_show)) - min(min(O_show)));
            imwrite(O_show,sprintf('./dbg_CS_PM3D/iter_%03d.png',ss));
%             figure(1),imshow(abs(u0_est),[]), title('CSGI');
        end %% END of MAIN ALGORITHM, ss
        O = u0_est;
        
end

% if params.norm % normalization
%     O = (O - min(min(O)))/(max(max(O)) - min(min(O)));
% end

end