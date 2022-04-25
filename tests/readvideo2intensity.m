function [ I ] = readvideo2intensity(videopath)
%READVIDEO2INTENSITY Read video frame-by-frame to get the intensity of the
%correlated measurement.
vid = VideoReader(videopath);
m = 1;
while hasFrame(vid)
    frame = readFrame(vid);
    if ndims(frame)==1 % monochrome
        I(m) = sum(frame(:));
    else % color
        [~,~,nch] = size(frame);
        for ich = 1:nch
            I(m,ich) = sum(frame(:,:,ich),'all'); % color channel order: RGB
        end
    m = m+1;
end

end