function WSC = LocalNormalConstraintMatrix(nMap, h, w, sig_n)
% Function for local normal(shading) constraint
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
    
    nMap_vec = reshape(nMap,h*w,3);
    nMap_p = nMap_vec( edge_sorted(:,1) , : );
    nMap_q = nMap_vec( edge_sorted(:,2) , : );
    dist = 1 - sum( nMap_p .* nMap_q , 2 );
    weight = exp(-dist.^2 / sig_n ^2);
    
    % WSC : h*w x h*w
    WSC = sparse(edge_sorted(:,1),edge_sorted(:,1),weight(:),h*w,h*w)...
        - sparse(edge_sorted(:,1),edge_sorted(:,2),weight(:),h*w,h*w)...
        + sparse(edge_sorted(:,2),edge_sorted(:,2),weight(:),h*w,h*w)...
        - sparse(edge_sorted(:,2),edge_sorted(:,1),weight(:),h*w,h*w); 
end