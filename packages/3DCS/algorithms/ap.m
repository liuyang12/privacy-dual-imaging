function [ xopt ] = ap( A,y,opts )
%AP Alternating projection (AP) method for single-pixel imaging
%   xopt=AP(A,y,opts) returns the optimal alternating projection result,
%   where A is m-by-n measurement matrix or sensing matrix, y is m-by-1
%   measurements, and opts are the parameters used in this solver.
%   [model]
%   
%   [default]
%     opts.maxiter                                        100
%   [refenrences] This implementation follows the Biomedical Optics Express
%   (BOE) paper by K. Guo in 2016.
%     Guo, K., Jiang, S., & Zheng, G. (2016). Multilayer fluorescence 
%       imaging on a single-pixel detector. Biomedical Optics Express, 
%       7(7), 2425-2431.
%   
%   Copyright(C) 2016-2017 by <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>.
%   Last modified Feb 9, 2017.
%   See also CS.
if nargin<3 || ~isfield(opts,'maxiter')
    opts.maxiter = 100; % [default] maximum iterations
end
maxiter = opts.maxiter; % maximum iterations

[m,n] = size(A);
% x = ones(1,n); % initial guess of x
x = zeros(1,n); % initial guess of x
% iteration
for it = 1:maxiter
    % each measurement
    for j = 1:m
        x1 = x.*A(j,:);
        x0 = x1;
        x1 = x1 - mean(x1) + y(j)/n;
        x = abs(x + A(j,:)./(max(A(j,:)).^2).*(x1-x0));
    end
end

xopt = x';
fprintf('AP ends at its %4d iteration.\n',it);

end

