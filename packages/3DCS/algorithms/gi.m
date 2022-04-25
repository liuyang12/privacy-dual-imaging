function [ giout ] = gi( sensmat,meas,params )
%GI Correlation-based algorithms in solving ghost imaging reconstruction
%problem
%   giout=GI(sensmat,meas,params) returns the N-by-1 vector of the sparse 
%   reconstruction result, where sensmat is the M-by-N sensing matrix or 
%   measurement matrix, M is the number of measurements and N is the 
%   Nyquist number of the orignal signal, meas is the M-by-1 measurement 
%   vector that satisies the compressiver sampling equation 
%   meas=sensmat*sig_out+nois, and nois here is the measurement noise on 
%   meas, params is a structure specifying the parameters of the solver, 
%   including the algorithm (method) with the corresponding options.
%   See also CS.
% 
%   Copyright(C) 2016-2017 by <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>.
%   Last modified Feb 9, 2017.

% [1] basic parameters
opts = []; % options for the reconstruction algorithms

% [2] apply correlation-based algorithms
switch lower(params.gimethod) % [params] gimethod
    case 'tgi'     % [2.1] traditional ghost imaging
        [M,N] = size(sensmat);
        measmean = 0;
        pattmean = zeros(1,N);
        for i = 1:M
            measmean = (measmean*(i-1)+meas(i))/i;
            pattmean = (pattmean*(i-1)+double(sensmat(i,:)))/i;
            obj = (obj*(i-1)+(meas(i)-measmean)*(double(sensmat(i,:))-pattmean))/i;
        end
        giout = obj';
        % [end] TGI
        %
        %
end

