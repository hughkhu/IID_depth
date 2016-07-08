function [W] = lleWeight(X,neighborhood)
% LLE ALGORITHM - ONLY Compute weights
% Adapted by HughKhu   2016-05-11
% [W] = lle(X,neighborhood)
% X = data as D x N matrix (D = dimensionality, N = #points)
% neighborhood = data as k x N (k = number of neighbor points)
% W = data as k x N
    [D,N] = size(X);
    fprintf(1,'LLE running on %d points in %d dimensions\n',N,D);
    K=size(neighborhood,1);%K nearest
    % SOLVE FOR RECONSTRUCTION WEIGHTS
    fprintf(1,'-->Solving for reconstruction weights.\n');
    if(K>D) 
      fprintf(1,'   [note: K>D; regularization will be used]\n'); 
      tol=1e-3; % regularlizer in case constrained fits are ill conditioned
    else
      tol=0;
    end
    
    W = zeros(K,N);
    for ii=1:N
       z = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K); % shift ith pt to origin
       C = z'*z;                                        % local covariance
       C = C + eye(K,K)*tol*trace(C);                   % regularlization (K>D)
       W(:,ii) = C\ones(K,1);                           % solve Cw=1
       W(:,ii) = W(:,ii)/sum(W(:,ii));                  % enforce sum(w)=1
    end
    fprintf(1,'Done.\n');
end
