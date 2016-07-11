addpath('filtering\\RGF');      % rolling guidance filter
% Read Data & Preprocessing (Structure-texture separation)
I = im2double(imread('images\\image0008.png'));
depth = double(imread('images\\raw_depth0008.png'))/1000.0;
%raw_depth in mm, here coverts to meter

% Structure-Texture Separation (L.Karacan,E. Erdem and A. Erdem.
% Structure Preserving Image Smoothing via Region Covariances. ACM
% Transactions on Graphics (Proceedings of SIGGRAPH Asia 2013), 32(6),
% November 2013) 
S = RollingGuidanceFilter(I, 3, 0.1, 4);
%%
slic = 1;
tic
if slic
    [reflectance, shading] = IID_slic(I, S, depth, 1);
else
    [reflectance, shading] = IID_adapted(I, S, depth, 1);
end
toc
figure;imshow(reflectance);
figure;imshow(shading);
