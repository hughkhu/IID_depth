addpath('filtering\\RGF');   % rolling guidance filter
%option: select IID_slic or IID_adapted
slic = 1;% set 0 if adapted.
%% path 
Path.DataSet='D:\NYU_RGBDdataset\';
Path.depFile='NYUdepths\';
Path.rgbFile='NYUimages\';
Path.imageFore='image';
Path.depthFore='raw_depth';
%id=1;
Path.rgbExt='.png';
if slic
    Path.wrFile = 'result-slic\';
else
    Path.wrFile = 'result-my_texture_remove\'; %my adaption
end
Path.ref='_reflectance';
Path.shd='_shading';
%% compute and output
% for id=[171,1,2,3,4,5,6,7]
for id=[172]
    tic;
    disp(['IID for ',num2str(id)]);
    depth_path=strcat(Path.DataSet,Path.depFile,...
                      Path.depthFore,sprintf('%04d',id),Path.rgbExt);
    depth=double(imread(depth_path))/1000;%in meter double

    rgb_path=strcat(Path.DataSet,Path.rgbFile,...
                      Path.imageFore,sprintf('%04d',id),Path.rgbExt);
    I=im2double(imread(rgb_path));
    S = RollingGuidanceFilter(I, 3, 0.1, 4);%remove texture
    %ÇÐ³ý°×±ß¿ò,crop_option = 1; compute
    if slic
        [reflectance, shading] = IID_slic(I, S, depth, 1);
    else
        [reflectance, shading] = IID_adapted(I, S, depth, 1);
    end
%  figure;imshow(reflectance);
%  figure;imshow(shading);
    disp(['Done: IID for ',num2str(id)]);
    %write
    shd_path=strcat(Path.DataSet,Path.wrFile,...
        sprintf('%04d',id),Path.shd,Path.rgbExt);
    ref_path=strcat(Path.DataSet,Path.wrFile,...
       sprintf('%04d',id),Path.ref,Path.rgbExt);
    imwrite(shading,shd_path);
    imwrite(reflectance,ref_path);

    tmptime=toc;
%     fid=fopen('_SLICruntimeNYU.txt','a+');
%     fprintf(fid,strcat(sprintf('%04d',id),': %4.2d seconds\t\n'),tmptime);
%     fclose(fid);
end
