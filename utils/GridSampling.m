function [ spl_pos, spl_nMap, spl_Points ] = GridSampling( nMap, Points, step)
% sample the image by square
% select one point whose variance of Normal is smallest in each square.
% written by HughKhu
% 2016-07-08
    [h,w,channel]=size(nMap);
    vMap=zeros(1,h*w);
    var_pad=[2,2];var_patch=[5,5];
    for i=1:channel
        vMap = vMap + var(im2col(padarray(nMap(:,:,i),var_pad,'symmetric'),var_patch,'sliding'));
    end
    vMap=reshape(vMap,[h,w]);
    nMap_vec=reshape(nMap,h*w,[]);
    Points_vec=reshape(Points,h*w,[]);
    %sampling，取normals，记录pos  %find(A=min(min(A)),1)
    %step=12;%使得整除
    spl_h=floor(h/step);
    spl_w=floor(w/step);
    spl_num=spl_h*spl_w;
    spl_pos=zeros(spl_num,1);%pos_r+pos_c*h
    spl_nMap=zeros(spl_num,3);%global-normal
    spl_Points=zeros(spl_num,3);%points
    count=0;
    for j=1:step:spl_w*step
        for i=1:step:spl_h*step
            count=count+1;
            winmat=vMap(i:i+step-1,j:j+step-1);
            [dr dc]=find(winmat==min(min(winmat)),1); % row and column of the minimum in the winmat patch
            spl_pos(count)=(j-1+dc-1)*h+(i+dr-1); % position in the whole image
            spl_nMap(count,:)=nMap_vec(spl_pos(count),:);
            spl_Points(count,:)=Points_vec(spl_pos(count),:);
        end
    end
end