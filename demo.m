addpath('filtering\\RGF');                     % rolling guidance filter
% Read Data & Preprocessing (Structure-texture separation)
I = im2double(imread('images\\image0171.png'));
depth = double(imread('images\\raw_depth0171.png'))/1000.0;
%raw_depth in mm, here coverts to meter

% Structure-Texture Separation (L.Karacan,E. Erdem and A. Erdem.
% Structure Preserving Image Smoothing via Region Covariances. ACM
% Transactions on Graphics (Proceedings of SIGGRAPH Asia 2013), 32(6),
% November 2013) 
S = RollingGuidanceFilter(I, 3, 0.1, 4);
%%
%{
tic
 [reflectance2, shading2] = IID_adapted(I, S, depth);
toc
 [reflectance, shading] = IID_slic(I, S, depth);
%}
[reflectance, shading] = myiidSLIC_noFilter(I, depth, 0.0001, 0.8, 0.5);
figure;imshow(reflectance);
figure;imshow(shading);
%{
imwrite(shading,'8s.png');
imwrite(reflectance,'8r.png');
%}