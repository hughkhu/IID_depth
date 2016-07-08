function WSC = LocalNormalConstraintMatrix(nMap, h, w, sig_n)
% function for local normal(shading) constraint
% written by HughKhu
% 2016-07-08
%% shangxia
    nMap_pt = nMap(2:h,:,:);
    nMap_qt = nMap(1:h-1,:,:);
    dist = 1 - sum (nMap_pt .* nMap_qt , 3);
    weight = exp(-dist.^2 / sig_n.^2);
    weight_t = reshape(weight,[],1);
    pos_pre = 1:h*w;
    pos_pre = reshape(pos_pre,h,w);
    pos_t = pos_pre-1;
    pos_pre = reshape(pos_pre(2:h,:),[],1);
    pos_t = reshape(pos_t(2:h,:),[],1);
    edge_t = [pos_pre,pos_t];
%% zuoyou
    nMap_pl = nMap(:,2:w,:);
    nMap_ql = nMap(:,1:w-1,:);
    dist = 1 - sum (nMap_pl .* nMap_ql , 3);
    weight = exp(-dist.^2 / sig_n.^2);
    weight_l = reshape(weight,[],1);
    pos_pre = 1:h*w;
    pos_pre = reshape(pos_pre,h,w);
    pos_l = pos_pre-h;
    pos_pre = reshape(pos_pre(:,2:w),[],1);
    pos_l = reshape(pos_l(:,2:w),[],1);
    edge_l = [pos_pre,pos_l];
%% zuoshang-youxia
    nMap_pdiag = nMap(2:h,2:w,:);
    nMap_qdiag = nMap(1:h-1,1:w-1,:);
    dist = 1 - sum (nMap_pdiag .* nMap_qdiag , 3);
    weight = exp(-dist.^2 / sig_n.^2);
    weight_diag = reshape(weight,[],1);
    pos_pre = 1:h*w;
    pos_pre = reshape(pos_pre,h,w);
    pos_diag = pos_pre-h-1;
    pos_pre = reshape(pos_pre(2:h,2:w),[],1);
    pos_diag = reshape(pos_diag(2:h,2:w),[],1);
    edge_diag = [pos_pre,pos_diag];
%% youshang-zuoxia
    nMap_pndiag = nMap(2:h,1:w-1,:);
    nMap_qndiag = nMap(1:h-1,2:w,:);
    dist = 1 - sum (nMap_pndiag .* nMap_qndiag , 3);
    weight = exp(-dist.^2 / sig_n.^2);
    weight_ndiag = reshape(weight,[],1);
    pos_pre = 1:h*w;
    pos_pre = reshape(pos_pre,h,w);
    pos_ndiag = pos_pre+h-1;
    pos_pre = reshape(pos_pre(2:h,1:w-1),[],1);
    pos_ndiag = reshape(pos_ndiag(2:h,1:w-1),[],1);
    edge_ndiag = [pos_pre,pos_ndiag];
%% sum	
    edge = [edge_t;edge_l;edge_diag;edge_ndiag];
    weight = [weight_t;weight_l;weight_diag;weight_ndiag];
    % sorting(optional)
    edge_weight = [min(edge(:,1),edge(:,2)),max(edge(:,1),edge(:,2)),weight];
    edge_weight = sortrows(edge_weight);
    % WSC : h*w x h*w
    WSC = sparse(edge_weight(:,1),edge_weight(:,1),edge_weight(:,3),h*w,h*w)...
        - sparse(edge_weight(:,1),edge_weight(:,2),edge_weight(:,3),h*w,h*w)...
        + sparse(edge_weight(:,2),edge_weight(:,2),edge_weight(:,3),h*w,h*w)...
        - sparse(edge_weight(:,2),edge_weight(:,1),edge_weight(:,3),h*w,h*w);    
end