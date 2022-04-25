function [para] = defparaconf(para)
%DEFPARACONF Default parameter configuration for WNNM-based image
%denoising.
%   para = DEFPARACONF(para) returns the parameters para by adding the
%   fields defined in the default parameters but undefined in the original
%   para.
%   See also PARACONFIG.
range = 255; % default range of the input signal for denoising
if isfield(para,'range'), range = para.range; end
dp = paraconfig(para.sigma*255/range); % default parameters
if isfield(para,'range')
    dp.nsigma = para.sigma/range;
else
    dp.nsigma = para.sigma;
end
if isfield(para,'patchsize') && ~isfield(para,'patchstep')
    para.patchstep = floor(para.patchsize/2-1);
end
fnames = fieldnames(dp);
for in = 1:length(fnames)
    fn = fnames{in};
    if ~isfield(para,fn)
        eval(sprintf('para.%s=dp.%s;',fn,fn));
    end
end

end

