function [ fringe ] = bin2grayfringe( fh,fv,phi,h,v,binning,opts )
%BIN2GRAYFRINGE Transfer binary fringe to grayscale fringe using digital
%dithering/halftoning or Fourier filtering.
%   
%   See also SINGLEPIXELQPISIM.

[X,Y] = meshgrid(1:h,1:v);
hc = floor(h/2)+1; % horizontal zero-frequency point after fftshift
vc = floor(v/2)+1; % vertical zero-frequency point after fftshift
P = 1/2+1/2*cos(2*pi/h*fh*X+2*pi/v*fv*Y+phi); % desired analog gray-scale fringe
if binning>1
    P = imresize(P,binning,'nearest');
end
switch opts.method
    case {'dither','halftone'} % digital dithering/halftoning method
        fringe = dither(P); % image dithering by applying Floyd-Steinberg's error diffusion algorithm.
    case 'filter' % Fourier filtering method
        if ~isfield(opts,'r_filt') % [default] pixel radius of the Fourier filter
            opts.r_filt = 1;
        end
        if ~isfield(opts,'cent_includ') % [default] center frequency included
            opts.cent_includ = 1;
        end
        if ~isfield(opts,'bin_dither') % [default] convert from grayscale to binary using dithering or not. 1 - dithering, 0 - not dithering.
            opts.bin_dither = 0;
        end
        r = opts.r_filt;
        if opts.bin_dither % using dithering
            Pbin = dither(P); % binary display on the digital micromirror device (DMD)
        else % not using dithering
            Pbin = double(P>0.5); % binary display on the digital micromirror device (DMD)
        end
        FP = fftshift(fft2(Pbin)); % fast discrete Fourier transform 
        if opts.cent_includ % center frequency included
            Ffilter = ((X-hc-fh).^2+(Y-vc-fv).^2<r^2) + ((X-hc).^2+(Y-vc).^2<r^2) + ((X-hc+fh).^2+(Y-vc+fv).^2<r^2) + ...
                      ((X-mod(hc+fh,h)).^2+(Y-mod(vc+fv,v)).^2<r^2) + ((X-mod(hc-fh,h)).^2+(Y-mod(vc-fv,v)).^2<r^2) + ...
                      ((X-hc-fh-h).^2+(Y-vc-fv-v).^2<r^2) + ((X-hc+fh-h).^2+(Y-vc+fv-v).^2<r^2) > 0.5; % initialization of the Fourier filter
        else % center frequency not included
            Ffilter = ((X-hc-fh).^2+(Y-vc-fv).^2<r^2) + ((X-hc+fh).^2+(Y-vc+fv).^2<r^2) + ...
                      ((X-mod(hc+fh,h)).^2+(Y-mod(vc+fv,v)).^2<r^2) + ((X-mod(hc-fh,h)).^2+(Y-mod(vc-fv,v)).^2<r^2) + ...
                      ((X-hc-fh-h).^2+(Y-vc-fv-v).^2<r^2) + ((X-hc+fh-h).^2+(Y-vc+fv-v).^2<r^2) > 0.5; % initialization of the Fourier filter
        end
        FPfilt = FP.*double(Ffilter);
        fringe = ifft2(ifftshift(FPfilt));
    otherwise
        error('Unsupported binary to grayscale fringe method %s.',opts.method);
end

end

