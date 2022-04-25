function [ im ] = imreadallfmt( impath )
%IMREADALLFMT read images of all formats
%   im=IMREADALLFMT(impath) returns the image matrix under such impath.
%   See also IMREAD.
imextensions = {'','.bmp','.jpg','.tif','.tiff','.cur','.gif','.hdf4','.ico','.pbm','.pcx','.pgm','.png','.ppm','.ras','.xwd'};
for i = 1:length(imextensions)
    impathfmt = [impath imextensions{i}];
    if exist(impathfmt,'file')
        if i == length(imextensions)
            error('specify image format and use imread instead.\n');
        else
            break;
        end
    end
end
im = imread(impathfmt); % read image

end

