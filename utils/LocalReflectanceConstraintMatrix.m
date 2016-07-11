function [WRC,consVecCont]=LocalReflectanceConstraintMatrix(chrom,S,h,w,sig_c,sig_i)
% Function for local reflectance constraint
% Author: HughKhu
% Version: 2.0
% Date: 2016-07-11
    pos_ori = reshape( 1 : h*w, h, w );
    % top-bottom
    pos_t1 = reshape(pos_ori(1:h-1,:),[],1);
    pos_t2 = reshape(pos_ori(2:h,:),[],1);
    edge_t = [pos_t1,pos_t2];
    % left-right
    pos_l1 = reshape(pos_ori(:,1:w-1),[],1);
    pos_l2 = reshape(pos_ori(:,2:w),[],1);
    edge_l = [pos_l1,pos_l2];
    % diagonal
    pos_diag1 = reshape(pos_ori(1:h-1,1:w-1),[],1);
    pos_diag2 = reshape(pos_ori(2:h,2:w),[],1);
    edge_diag = [pos_diag1,pos_diag2];
    % back-diagonal
    pos_bdiag1 = reshape(pos_ori(2:h,1:w-1),[],1);
    pos_bdiag2 = reshape(pos_ori(1:h-1,2:w),[],1);
    edge_bdiag = [pos_bdiag1,pos_bdiag2];
    
    edge = [edge_t; edge_l; edge_diag; edge_bdiag];
    edge_sorted = sortrows(edge);% sorting(optional)
    
    chrom_vec = reshape(chrom,h*w,3);
    chrom_p = chrom_vec( edge_sorted(:,1) , : );
    chrom_q = chrom_vec( edge_sorted(:,2) , : );
    dist = 2.0 * ( 1 - sum( chrom_p .* chrom_q , 2 ) );
    
    S_vec = reshape(S,h*w,3);
    S_p = S_vec( edge_sorted(:,1) , : );
    S_q = S_vec( edge_sorted(:,2) , : );
    S2p = max( sqrt( sum( S_p .^ 2 , 2 ) ) , 0.0001 );
    S2q = max( sqrt( sum( S_q .^ 2 , 2 ) ) , 0.0001 );
    
    weight = 1 + exp( - S2p .^ 2 / sig_i ^2 - S2q .^ 2 / sig_i ^2);
    weight = weight .* exp( -dist .^2 / sig_c ^2);
    
    dS = log( S2p ) - log( S2q );
    dweight = dS .* weight;
    % edge_weight=[edge_sorted,weight,dweight];
    
    % WSC : h*w x h*w
    WRC = sparse(edge_sorted(:,1),edge_sorted(:,1),weight(:),h*w,h*w)...
        - sparse(edge_sorted(:,1),edge_sorted(:,2),weight(:),h*w,h*w)...
        + sparse(edge_sorted(:,2),edge_sorted(:,2),weight(:),h*w,h*w)...
        - sparse(edge_sorted(:,2),edge_sorted(:,1),weight(:),h*w,h*w); 
     
    npairs=size(edge,1);
    tmpmap = sparse(edge_sorted(:,1),1:npairs,ones(npairs,1),h*w,npairs);
    consVecCont = tmpmap * dweight;
    tmpmap = sparse(edge_sorted(:,2),1:npairs,ones(npairs,1),h*w,npairs);
    consVecCont = consVecCont - tmpmap * dweight;
end