function [ res_r, res_s ] = IID_adapted( I, S, depth, crop_option )
%   It require RGB image, its texture-removed image, and aligned depth image as inputs.
%   It decompose input image into reflectance and shading image.
%   Parameters
%   I: input RGB image (double)
%   S: texture-removed RGB image
%   depth: Kinect depth image (in meter).
%   crop_option: crop the image or not. 1 for yes, 0 for no.
%   ==========
%   The Code is created adapted from the method described in the following paper:
%   [1] "Intrinsic Image Decomposition Using Structure-Texture Separation
%   and Surface Normals", Junho Jeon, Sunghyun Cho, Xin Tong, Seungyong
%   Lee, European Conference on Computer Vision (ECCV), 2014

    disp('Previous Processing...');
    run('D:\MatlabWorks\vlfeat-0.9.19/toolbox/vl_setup.m');% for KD-Tree 
    addpath('utils'); % basic utilities
    %some para,not recommend to change.
    sigma_c = 0.0001;
    sigma_i = 0.8;
    sigma_n = 0.5;
  
    D = depth;
    [h, w, ~] = size(I);
    D = smooth_d(S, double(D), ones(h, w)); % RGB-D Joint bilateral filtering (Code from Qifen Chen and Vladlen Koltun, ICCV 2013)
    Points=getVectors(size(D,1),size(D,2));
	Points=Points .* D(:,:,[1 1 1]);
    % Normal Map Estimation
    [nx, ny, nz]=surfnorm(Points(:,:,1),Points(:,:,2),Points(:,:,3));
    nMap=cat(3,nx,ny,nz);
    % Cropping (Kinect RGB image has a white padding)
    if crop_option
        S = S(7:h-6, 9:w-8, :);
        I = I(7:h-6, 9:w-8, :);
        Points = Points(7:h-6, 9:w-8, :);     
    end
    [Points, nMap] = DenoisePoints(Points, nMap); % Bilateral Mesh Denoising
    [h, w, ~] = size(I);
    N = h*w;
    % Minimum Patch Normal Variance Sub-Sampling
    [ sample_pos, sample_nMap, sample_Points ] = GridSampling( nMap, Points, 12);
%% Sub-Sampled non-local LLE Constraint
    k=50;%num of neighbors
    feature=sample_nMap';%feature vector (n_x,n_y,n_z)' 3*sample_num
    idx=double(vl_kdtreequery(vl_kdtreebuild(feature),feature,feature,'numneighbors',k+1));
    neig_idx=(idx(2:k+1,:))';%sample_num*k
    neig_weight=(lleWeight(sample_nMap',neig_idx'))';%LLE

    edge=[ repmat(sample_pos,k,1),sample_pos(neig_idx(:)) ];
    NonLocalLLEMatrix = sparse(edge(:,1),edge(:,2),reshape(neig_weight,[],1),N,N);
    LLENORMAL = NonLocalLLEMatrix;
%% Sub-Sampled local LLE Constraint  
    k=50;%num of neighbors
    feature=sample_nMap';%feature vector (x,y,z,n_x,n_y,n_z)' 6*sample_num
    feature(4:6,:) = sample_Points';
    idx2=double(vl_kdtreequery(vl_kdtreebuild(feature),feature,feature,'numneighbors',k+1));
    neig_idx2=(idx2(2:k+1,:))';%sample_num*k
    neig_weight2=(lleWeight(sample_nMap',neig_idx2'))';%LLE

    edge=[ repmat(sample_pos,k,1), sample_pos(neig_idx2(:)) ];
    LocalLLEMatrix = sparse(edge(:,1),edge(:,2),reshape(neig_weight2,[],1),N,N);
    LLEGRID = LocalLLEMatrix;   
%% Matting Laplacian matrix 
    % Compute propagation weights (matting Laplacian)
    disp('computing laplacian matrix.');
    L_S = getLaplacian1(S, zeros(h, w), 0.1^5);    
%% Local Retinex Constraints   
    C = getChrom(S);
    % Compute local normal(shading) constraint (continuous similarity weight)
    WSC=LocalNormalConstraintMatrix(nMap,h,w,sigma_n);
    % Compute local reflectance constraint (continuous similarity weight)
    [WRC,consVecCont]=LocalReflectanceConstraintMatrix(C,S,h,w,sigma_c,sigma_i);
%% Optimization
    spI = speye(N, N);    
    mk = zeros(N, 1);
    mk(sample_pos) = 1;
    mask = spdiags(mk, 0, N, N); % Subsampling mask for local LLE
    mk = zeros(N, 1);
    mk(sample_pos) = 1;
    mask2 = spdiags(mk, 0, N, N); % Subsampling mask for non-local LLE
    
    A = 4 * WRC + 1 * mask * (spI - LLEGRID) + 1 * mask2 * (spI - LLENORMAL) + 0.025 * WSC + L_S;
%   A = 4 * WRC + 1 * mask * (spI - LLEGRID) + 1 * mask2 * (spI - LLENORMAL) + 0.025 * WSC;
    b = 4 * consVecCont;    

    disp('Optimizing the system...');
    newS = pcg(A, b, 1e-3, 10000, [], []);
    % Visualization and Saving Results
    res_s = reshape(exp(newS), [h w])/2;
    res_r = I ./ repmat(res_s, [1 1 3]) /2;
end


