function [WRC,consVecCont]=LocalReflectanceConstraintMatrix(chrom,S,h,w,sig_c,sig_i)
% function for local reflectance constraint
% written by HughKhu
% 2016-07-08
%% shangxia	 
    S_pt = S(2:h,:,:);
    S_qt = S(1:h-1,:,:);
    chrom_pt = chrom(2:h,:,:);
    chrom_qt = chrom(1:h-1,:,:);
    S2p = max(sqrt(sum(S_pt .^ 2 , 3)),0.0001);
    S2q = max(sqrt(sum(S_qt .^ 2 , 3)),0.0001);
    dist=2.0*(1-(sum(chrom_pt.*chrom_qt,3)));
    weight = 1 + exp(-S2p.^2/sig_i^2 - S2q.^2/ sig_i^2);
    weight = weight .* exp(-dist.^2/sig_c.^2);
    S2p = log(S2p);
    S2q = log(S2q);
    dS = S2p - S2q;
    dweight = dS .* weight;
    weight_t = reshape(weight,[],1);
    dweight_t = reshape(dweight,[],1);
    pos_pre = 1 : h*w;
    pos_pre = reshape(pos_pre,h,w);
    pos_t = pos_pre - 1;
    pos_pre = reshape(pos_pre(2:h,:),[],1);
    pos_t = reshape(pos_t(2:h,:),[],1);
    edge_t = [pos_pre,pos_t];
%% zuoyou
    S_pl = S(:,2:w,:);
    S_ql = S(:,1:w-1,:);
    chrom_pl = chrom(:,2:w,:);
    chrom_ql = chrom(:,1:w-1,:);
    S2p = max(sqrt(sum(S_pl .^ 2 , 3)),0.0001);
    S2q = max(sqrt(sum(S_ql .^ 2 , 3)),0.0001);
    dist = 2.0*(1-(sum(chrom_pl.*chrom_ql,3)));
    weight = 1 + exp(-S2p.^2/sig_i^2 - S2q.^2/ sig_i^2);
    weight = weight .* exp(-dist.^2/sig_c.^2);
    S2p = log(S2p);
    S2q = log(S2q);
    dS = S2p - S2q;
    dweight = dS .* weight;
    weight_l = reshape(weight,[],1);
    dweight_l = reshape(dweight,[],1);
    pos_pre = 1 : h*w;
    pos_pre = reshape(pos_pre,h,w);
    pos_l = pos_pre - h;
    pos_pre = reshape(pos_pre(:,2:w),[],1);
    pos_l = reshape(pos_l(:,2:w),[],1);
    edge_l = [pos_pre,pos_l];
%% zuoshang-youxia
    S_pdiag = S(2:h,2:w,:);
    S_qdiag = S(1:h-1,1:w-1,:);
    chrom_pdiag = chrom(2:h,2:w,:);
    chrom_qdiag = chrom(1:h-1,1:w-1,:);
    S2p = max(sqrt(sum(S_pdiag .^ 2 , 3)),0.0001);
    S2q = max(sqrt(sum(S_qdiag .^ 2 , 3)),0.0001);
    dist = 2.0*(1-(sum(chrom_pdiag.*chrom_qdiag,3)));
    weight = 1 + exp(-S2p.^2/sig_i^2 - S2q.^2/ sig_i^2);
    weight = weight .* exp(-dist.^2/sig_c.^2);
    S2p = log(S2p);
    S2q = log(S2q);
    dS = S2p - S2q;
    dweight = dS .* weight;
    weight_diag = reshape(weight,[],1);
    dweight_diag = reshape(dweight,[],1);
    pos_pre = 1 : h*w;
    pos_pre = reshape(pos_pre,h,w);
    pos_diag = pos_pre - h - 1;
    pos_pre = reshape(pos_pre(2:h,2:w),[],1);
    pos_diag = reshape(pos_diag(2:h,2:w),[],1);
    edge_diag = [pos_pre,pos_diag];
%% youshang-zuoxia
    S_pndiag = S(2:h,1:w-1,:);
    S_qndiag = S(1:h-1,2:w,:);
    chrom_pndiag = chrom(2:h,1:w-1,:);
    chrom_qndiag = chrom(1:h-1,2:w,:);
    S2p = max(sqrt(sum(S_pndiag .^ 2 , 3)),0.0001);
    S2q = max(sqrt(sum(S_qndiag .^ 2 , 3)),0.0001);
    dist = 2.0*(1-(sum(chrom_pndiag.*chrom_qndiag,3)));
    weight = 1 + exp(-S2p.^2/sig_i^2 - S2q.^2/ sig_i^2);
    weight = weight .* exp(-dist.^2/sig_c.^2);
    S2p = log(S2p);
    S2q = log(S2q);
    dS = S2p - S2q;
    dweight = dS .* weight;
    weight_ndiag = reshape(weight,[],1);
    dweight_ndiag = reshape(dweight,[],1);
    pos_pre = 1 : h*w;
    pos_pre = reshape(pos_pre,h,w);
    pos_ndiag = pos_pre+h-1;
    pos_pre = reshape(pos_pre(2:h,1:w-1),[],1);
    pos_ndiag = reshape(pos_ndiag(2:h,1:w-1),[],1);
    edge_ndiag = [pos_pre,pos_ndiag];
%% sum
    edge=[edge_t;edge_l;edge_diag;edge_ndiag];
    weight=[weight_t;weight_l;weight_diag;weight_ndiag];
    dweight=[dweight_t;dweight_l;dweight_diag;dweight_ndiag];
    edge_weight=[edge,weight,dweight];
    % not sorting: 左小对第四列有正负错误影响
    % edge_weight=[min(edge(:,1),edge(:,2)),max(edge(:,1),edge(:,2)),weight,dweight];
    % edge_weight=sortrows(edge_weight);

    WRC = sparse(edge_weight(:,1),edge_weight(:,1),edge_weight(:,3),h*w,h*w)...
        - sparse(edge_weight(:,1),edge_weight(:,2),edge_weight(:,3),h*w,h*w)...
        + sparse(edge_weight(:,2),edge_weight(:,2),edge_weight(:,3),h*w,h*w)...
        - sparse(edge_weight(:,2),edge_weight(:,1),edge_weight(:,3),h*w,h*w);   
    
    npairs=size(edge_weight,1);
    tmpmap = sparse(edge_weight(:,1),1:npairs,ones(npairs,1),h*w,npairs);
    consVecCont = tmpmap * dweight;
    tmpmap = sparse(edge_weight(:,2),1:npairs,ones(npairs,1),h*w,npairs);
    consVecCont = consVecCont - tmpmap * dweight;
end