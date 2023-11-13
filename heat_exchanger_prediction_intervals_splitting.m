
%Find rule with the highest error that was not removed
nonminMSE = MSE;
nonminMSE(n_c < 2*median(n_c)) = 0; %Remove clusters that have a low number of samples
[maxsize, ic] = max(nonminMSE); 

%Split the cluster with the highest error, 
% something might not be right with it
if  maxsize > 0 

    %Compute the SVD of the cluster to get the principle vector
    [~,S,V] = svd(S_c(:,:,ic)/n_c(ic));

    %Split the number of samples 
    n_pq = n_c(ic);
    n_p = floor(n_c(ic)/2); %Če naredim tako rata k_split 1 problem
    n_q = n_pq - n_p;

    %Setup the number of 
    E_p = ones(n_p, n_z);
    E_q = ones(n_q, n_z);
    E_pq = ones(n_pq, n_z);

    %Compute split mean
    mi_p = mi(ic,:) + k_split*sqrt(S(1))*V(:,1)';
    mi_q = mi(ic,:) - k_split*sqrt(S(1))*V(:,1)';

    %Setup splitting means matrixes for computations bellow
    M_q = diag(mi_q);
    M_p = diag(mi_p);
    M_pq = diag(mi(ic,:));

    %Compute the cluster covaraince matrix
    Sigma_p = (1/(n_q+n_p-2))*(((n_pq-1)/(n_pq))*S_c(:,:,ic) + ...
        M_pq'*(E_pq')*E_pq*M_pq - M_p'*(E_p')*E_p*M_p - M_q'*(E_q')*E_q*M_q);

    %Compute the nonnormalized cluster covariance matrix
    S_p = n_p*Sigma_p;
    S_q = n_q*Sigma_p;

    %Check for debugging
    if det(S_p)<0
        disp('Error spliting, negative variance')
    end

    %Add the 1. new cluster
    n_c(c+1) = n_p;
    mi(c+1,:) = mi_p;
    S_c(:,:,c+1) = S_p;
    theta(:,c+1) = theta(:,ic);
    MSE(c+1) = MSE(ic);
    b_fuzzy(c+1) = b_fuzzy(ic);
    P(:,:,c+1) = P(:,:,ic);

    %Compute  1. new cluster projection
%     Sigma(:,:,c+1) = S_c(:,:,c+1)/n_c(c+1);
%     [~,S,V] = svd(Sigma(:,:,c+1));
%     delta(:,:,c+1) = max(abs(V(1,:)'.*diag(sqrt(S))));
    phi_c(c+1,:) = phi_c(ic,:);

    %Add the 2. new cluster
    n_c(c+2) = n_q;
    mi(c+2,:) = mi_q;
    S_c(:,:,c+2) = S_q;
    theta(:,c+2) = theta(:,ic);
    MSE(c+2) = MSE(ic);
    b_fuzzy(c+2) = b_fuzzy(ic);
    P(:,:,c+2) = P(:,:,ic);
    phi_c(c+2,:) = phi_c(ic,:);

    %Compute  2. new cluster projection
%     Sigma(:,:,c+2) = S_c(:,:,c+2)/n_c(c+2);
%     [~,S,V] = svd(Sigma(:,:,c+2));
%     delta(:,:,c+2) = max(abs(V(1,:)'.*diag(sqrt(S))));

    %Remove the original cluster that was split
    n_c(ic) = [];
    mi(ic,:) = [];
    S_c(:,:,ic) = [];
    theta(:,ic) = [];
    MSE(ic) = [];
    b_fuzzy(ic) = [];
    P(:,:,ic) = [];
    Sigma(:,:,ic) = [];
    delta(:,:,ic) = [];


    %Increase the rule counter
    c = c+1;

end
