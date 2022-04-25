function [ im_norm ] = imnorm( im_orig )
%IMNORM Maximum-minimum normalization of the 2-d image to the range of [0,1].

if max(im_orig(:))==min(im_orig(:))
    im_norm = im_orig;
else
im_double = double(im_orig);
im_norm =    (im_double - min(im_double(:)))...
           / (max(im_double(:)) - min(im_double(:)));
end

end

