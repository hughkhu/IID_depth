function [ res_r, res_s ] = IID_slic( I, S, depth, crop_option)
%   Matlab implementation of My Undergraduate thesis about 
%   Intrinsic Image Decomposition using Depths.
%
%   Parameters
%   I: input RGB image (double)
%   S: texture-removed RGB image 
%   depth: Kinect depth image (in meter).
%   crop_option: crop the image or not. 1 for yes, 0 for no.
%   sigma_c, sigma_i, sigma_n: modification is not recommended.
%   ==========
%   The Code is created by reference to the following paper:
%   [1] "Intrinsic Image Decomposition Using Structure-Texture Separation
%   and Surface Normals", Junho Jeon, Sunghyun Cho, Xin Tong, Seungyong
%   Lee, European Conference on Computer Vision (ECCV), 2014
%
%   Author: HughKhu
%   Date  : 2016-07-11
%   Version : 1.0 
    disp('Previous Processing...');
    run('D:\MatlabWorks\vlfeat-0.9.19/toolbox/vl_setup.m');
    addpath('utils');            % basic utilities
    %some para,not recommend to change.
    sigma_c = 0.0001;
    sigma_i = 0.8;
    sigma_n = 0.5;
    
    D = depth;
    [h, w, ~] = size(I);
    D = smooth_d(S, double(D), ones(h, w)); % RGB-D Joint bilateral filtering (Code from Qifen Chen and Vladlen Koltun, ICCV 2013)
    Points = getVectors(size(D,1),size(D,2));
	Points = Points .* D(:,:,[1 1 1]);
    % Normal Map Estimation
    [nx, ny, nz] = surfnorm(Points(:,:,1),Points(:,:,2),Points(:,:,3));
    nMap = cat(3,nx,ny,nz);
    % Cropping (Kinect RGB image has a white padding)
    if crop_option
        S = S(7:h-6, 9:w-8, :);
        I = I(7:h-6, 9:w-8, :);
        Points = Points(7:h-6, 9:w-8, :);           % cropping the image
    end
    [Points, nMap] = DenoisePoints(Points, nMap); % Bilateral Mesh Denoising
    [h, w, ~] = size(I);
    N = h*w;
    
    C = getChrom(S);
    nMap_vec = reshape(nMap,N,[]);%N*3
    Points_vec = reshape(Points,N,[]);%N*3
%% Local Retinex Constraints
    disp('Local Retinex Constraints Computing...');
    % Compute local normal(shading) constraint (continuous similarity weight)
    WSC = LocalNormalConstraintMatrix(nMap,h,w,sigma_n);
    % Compute local reflectance constraint (continuous similarity weight)
    [WRC,consVecCont] = LocalReflectanceConstraintMatrix(C,S,h,w,sigma_c,sigma_i);

%% SLIC sample
    disp('SLIC Grouping...');
    %[labels, numlabels] = slicmex(I,2000,20);%numlabels is the same as number of superpixels
    [labels, numlabels] = slicmex(im2uint8(I),2000,20);%numlabels is the same as number of superpixels
    labels = labels+1;
    % figure,imagesc(labels);
    % colormap jet;
    sample_num = numlabels;
    sample_pos = zeros(sample_num,1);%pos_r+pos_c*h
    sample_nMap = zeros(sample_num,3);%global-normal
    sample_Points = zeros(sample_num,3);%points
    edgeSLIC = [];%inner group edgs
    for i = 1:numlabels
        labelcurr = find(labels==i); %a position column 
        pnum = numel(labelcurr);
        % sample and record
        sample_pos(i) = labelcurr(randi(pnum,1,1));
        sample_nMap(i,:) = nMap_vec(sample_pos(i),:);
        sample_Points(i,:) = Points_vec(sample_pos(i),:);
        % inner group edges
        edge = [repmat(labelcurr,10,1),labelcurr(randi(pnum,10*pnum,1))];
        edge = edge(find(edge(:,1)~=edge(:,2)),:);%去除相等的行
        edge = unique([min(edge(:,1),edge(:,2)) max(edge(:,1),edge(:,2))],'rows');%去掉重复的行
        edgeSLIC = [edgeSLIC;edge];% to be revised
    end
%% super pixel inner Group chrom
    disp('SLIC inner Group Constraint Computing...');   
    Chrom_vec = reshape(C,N,3);
    S_vec = reshape(S,N,3);
    dist = 1 - sum( Chrom_vec(edgeSLIC(:,1),:) .* Chrom_vec(edgeSLIC(:,2),:) , 2 );
    S2p = max ( sum( S_vec(edgeSLIC(:,1)).^2 , 2) , 1e-8);
    S2q = max ( sum( S_vec(edgeSLIC(:,2)).^2 , 2) , 1e-8);
    Cweight = 1 + exp(-(S2p+S2q) / sigma_i^2);
	Cweight = Cweight .* exp(-(2 * dist).^2 / sigma_c.^2);
    dS = log(sqrt(S2p)) - log(sqrt(S2q));
    dweight=dS.*Cweight;
    GroupChromMat = sparse(edgeSLIC(:,1),edgeSLIC(:,1),Cweight,N,N)...
                  - sparse(edgeSLIC(:,1),edgeSLIC(:,2),Cweight,N,N)...
                  + sparse(edgeSLIC(:,2),edgeSLIC(:,2),Cweight,N,N)...
                  - sparse(edgeSLIC(:,2),edgeSLIC(:,1),Cweight,N,N);
    GroupChrom_b =  sparse(edgeSLIC(:,1),1,dweight,N,1)...
                 -  sparse(edgeSLIC(:,2),1,dweight,N,1);      
%% Compute Sub-Sampled non-local LLE Constraint
    disp('non-local LLE Computing...');
    k=50;%num of neighbors
    feature=sample_nMap';%feature vector (n_x,n_y,n_z)' 3*sample_num
    idx=double(vl_kdtreequery(vl_kdtreebuild(feature),feature,feature,'numneighbors',k+1));
    neig_idx=(idx(2:k+1,:))';%sample_num*k
    neig_weight=(lleWeight(sample_nMap',neig_idx'))';%LLE
    
    edge=[reshape(repmat(sample_pos,1,k),[],1),sample_pos(neig_idx(:))];
    NonLocalLLEMatrix=sparse(edge(:,1),edge(:,2),reshape(neig_weight,[],1),N,N);
    LLENORMAL=NonLocalLLEMatrix;
%% Compute Sub-Sampled local LLE Constraint  
    disp('local LLE Computing...');
    k=50;%num of neighbors
    feature=sample_nMap';%feature vector (x,y,z,n_x,n_y,n_z)' 6*sample_num
    feature(4:6,:) = sample_Points';
    idx2=double(vl_kdtreequery(vl_kdtreebuild(feature),feature,feature,'numneighbors',k+1));
    neig_idx2=(idx2(2:k+1,:))';%sample_num*k
    neig_weight2=(lleWeight(sample_nMap',neig_idx2'))';%LLE
    
    edge=[reshape(repmat(sample_pos,1,k),[],1),sample_pos(neig_idx2(:))];
    LocalLLEMatrix=sparse(edge(:,1),edge(:,2),reshape(neig_weight2,[],1),N,N);
    LLEGRID=LocalLLEMatrix;   
%% Optimization
    spI = speye(N, N);    
    mk = zeros(N, 1);
    mk(sample_pos) = 1;
    mask = spdiags(mk, 0, N, N); % Subsampling mask for local LLE
    mk = zeros(N, 1);
    mk(sample_pos) = 1;
    mask2 = spdiags(mk, 0, N, N); % Subsampling mask for non-local LLE
    
    A = 4 * WRC + 1 * mask * (spI - LLEGRID) + 1 * mask2 * (spI - LLENORMAL) ...
      + 0.025 * WSC+ GroupChromMat;
    b = 4 * consVecCont + GroupChrom_b;    
%  tic
    disp('Optimizing the system...');
    newS = pcg(A, b, 1e-3, 10000, [], []);
    % Visualization and Saving Results
    res_s = reshape(exp(newS), [h w])/2;
    res_r = I ./ repmat(res_s, [1 1 3]) /2;
%  toc
end