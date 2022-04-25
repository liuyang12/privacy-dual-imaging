function imbin = imbinning ( im, binsize )
%IMBINNING Down-sampling an image using pixel binning, that is grouping the
%binning pixel(s) and summing (or averaging) them into super-pixels.
%  IMBIN = IMBINNING(IM, BINSIZE) returns the down-sampled image with pixel
%  binning, where BINSIZE is the binning size. 
%  
%  Note that this is like image (two-dimensional, or 2D) convolution with a
%  flat convolution kernel and the stride size same as the kernel size.
%
if length(binsize) < 2
    binsize = [binsize binsize]; % the same binning size for both dimensions
end

rowbin = binsize(1);
colbin = binsize(2);

imbin = zeros(size(im,1)/rowbin, size(im,2)/colbin, size(im,3));
for ir = 1:rowbin
    for ic = 1:colbin
        imbin = imbin + im(ir:rowbin:end,ic:colbin:end,:); % sum
    end
end

% imbin = imbin/(rowbin*colbin); % average

end

