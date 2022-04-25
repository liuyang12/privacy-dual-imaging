function [ im ] = readrawimage( im_path, opts )
%READRAWIMAGE Read .RAW image file
%   img=READRAWIMAGE(img_path,params) returns the image img in the 
%   destinated image path img_path formatted as .RAW.
%   [0] parameters
%       opts.rows         number of rows
%       opts.cols         number of columns
%       opts.format       format of the .RAW file
%   See also FREAD, IMREAD.

rows   = opts.rows; % number of rows
cols   = opts.cols; % number of columns
format = opts.format; % format of the .RAW image file

fin = fopen(im_path,'r'); % read the .RAW image file

I = fread(fin,[cols rows],format);
im = I';

fclose(fin);
end

