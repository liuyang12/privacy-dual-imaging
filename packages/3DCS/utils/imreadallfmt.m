function [ im ] = imreadallfmt( impath )
%IMREADALLFMT read images of all formats
%   im=IMREADALLFMT(impath) returns the image matrix under such impath.
%   See also IMREAD.
exts = {'.jpg', '.jpeg', '.png', '.bmp', '.tif', '.tiff'};

KNOWNIMFILE = false;
for i = 1:length(exts)
    fext = exts{i};
    impathfmt = [impath fext];
    if exist(impathfmt, 'file')
        KNOWNIMFILE = true;
        break;
    end
end

if KNOWNIMFILE
    im = imread(impathfmt); % read image
else
    error('specify image format and use imread instead.\n');
end

end

