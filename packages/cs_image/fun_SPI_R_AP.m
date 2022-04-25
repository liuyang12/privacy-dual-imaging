function [im_r] = fun_SPI_R_AP(patterns, measurements)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Liheng Bian, June 25, 2016
% Contact: lihengbian@gmail.com
% This function implements the singel pixel imaging reconstruction using the alternating projection method proposed in 
% "Kaikai Guo, Shaowei Jiang, and Guoan Zheng, Multilayer fluorescence imaging on a single-pixel detector, Biomedical Optics Express, 7, 7, 2425 (2016).".

% Inputs:
% patterns: illumination patterns (pixels * pixels * pattern numbers)
% measurements: single pixel measurements (vector)

% Outputs:
% im_r: reconstructed image  (pixels * pixels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im_r = ones(size(patterns,1),size(patterns,2));
for j = 1 : 5*30
    % begin iteration
    for i = 1 : size(patterns,3)
        temp_im = im_r .* patterns(:,:,i);
        temp_im1 = temp_im;    
        temp_im = temp_im - mean(mean(temp_im)) + measurements(i)/numel(temp_im);
        
        im_r = abs(im_r + (patterns(:,:,i))./((max(max(patterns(:,:,i)))).^2) .* (temp_im - temp_im1));
    end;
%     if mod(j,30) == 0
%         fprintf(['AP iteration ' num2str(j) '\n']);
        
% %         figure;
% %         imshow(abs(im_r),[],'InitialMagnification',1000);
% %         title(['Recovered im using AP method (iteration ' num2str(j) ')']);
% %         pause(0.1);
%     end
end
fprintf('AP iteration %.3d \n',j);
end