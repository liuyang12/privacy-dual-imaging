clear; clc; 
% close all;

H = 48; W = 48; % image size
h = 9; w = 9; % kernel size

imrec = [(1/4-1/24)*W, 1/4*H, (1/4+1/24)*W, 3/4*H;
         (1/2-1/24)*W, 1/4*H, (1/2+1/24)*W, 7/12*H;
         (3/4-1/24)*W, 1/4*H, (3/4+1/24)*W, 3/4*H];
     
kernrec = [3 3, 7 5;
           5 6, 5 9];

im = ones([H,W]);
kernel = zeros([h,w]);

for i = 1:size(imrec,1)
    p = round(imrec(i,1:2));
    q = round(imrec(i,3:4));
    
    im = make_rectangle(im, p, q, 0);
end

for k  =1:size(kernel,1)
    p = kernrec(k,1:2);
    q = kernrec(k,3:4);
    
    kernel = make_rectangle(kernel, p, q, 1);
end


function [im] = make_rectangle(im, p, q, color)
im(p(2):q(2),p(1):q(1)) = color;
end