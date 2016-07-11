addpath('filtering\\RGF');   % rolling guidance filter
%option: select IID_slic or IID_adapted
slic = 1;% set 0 if adapted.
%% path 
ALLfile={'alley_1','alley_2','ambush_2','ambush_4','ambush_5',...
         'ambush_6','ambush_7','bamboo_1','bamboo_2','bandage_1',...
         'bandage_2','cave_2','cave_4','market_2','market_5',...
         'market_6','mountain_1','shaman_2','shaman_3','sleeping_1',...
         'sleeping_2','temple_2','temple_3'};
         
Path.DataSet='D:\MPI-Sintel_Dataset\';
Path.depFile='depth\';
Path.rgbFile='clean\';

Path.frameFore='\frame_';
%id=1;
Path.depExt='.dpt';
Path.rgbExt='.png';
if slic
    Path.wrFile = 'result-slic\';
else
    Path.wrFile = 'result-my_texture_remove\'; %my adaption
end
% Path.wrFile = 'result-texture_remove\';%原代码需image的高、宽在sample时被整除，不适用于Sintel
Path.ref='_reflectance';
Path.shd='_shading';
%% compute and output
for j=23
    Path.clipFile=ALLfile{j};
    for id=12
        tic;
        disp(['IID for ',num2str(j),'-',num2str(id)]);
        depth_path=strcat(Path.DataSet,Path.depFile,Path.clipFile,...
                          Path.frameFore,sprintf('%04d',id),Path.depExt);
        depth=depth_read(depth_path);%in meter double

        rgb_path=strcat(Path.DataSet,Path.rgbFile,Path.clipFile,...
                          Path.frameFore,sprintf('%04d',id),Path.rgbExt);
        I=im2double(imread(rgb_path));

        S = RollingGuidanceFilter(I, 3, 0.1, 4);%remove texture
        %取消切除白边框, crop_option = 0; compute
        tic
        if slic
            [reflectance, shading] = IID_slic(I, S, depth, 0);
        else
            [reflectance, shading] = IID_adapted(I, S, depth, 0);
        end
        toc
        %{
            figure;imshow(reflectance);
            figure;imshow(shading);
        %}
        disp(['Done: IID for ',num2str(j),'-',num2str(id)]);
        %write
        shd_path=strcat(Path.DataSet,Path.wrFile,...
            Path.clipFile,sprintf('_%04d',id),Path.shd,Path.rgbExt);
        ref_path=strcat(Path.DataSet,Path.wrFile,...
            Path.clipFile,sprintf('_%04d',id),Path.ref,Path.rgbExt);
        imwrite(shading,shd_path);
        imwrite(reflectance,ref_path);

        tmptime=toc;
%         fid=fopen('D:\MatlabWorks\computeLMSE\_SLICruntime.txt','a+');
%         fprintf(fid,strcat(Path.clipFile,sprintf('_%04d',id),': %4.2d seconds\t\n'),tmptime);
%         fclose(fid);
    end
end