
%Compute overlapping
olap = Inf(c,c);
for i = 1:c
    %Compute overlapping measure: Bhattacharyya distance
    olap(i, :) = bhattacharyya(S_c, n_c, mi, c, i);
end

%Compute consequence dissimilarity
ar_pq = Inf(c,c);
Kappa_pq = Inf(c,c);

for i = 1:c
    for j = i+1:c

        %Transfer function of the current model
        B_p = theta(1:m, i);
        A_p = [1,theta(m+1:m+n, i)'];
        G_pu = (sum(B_p))/(sum(A_p));

        %Transfer function of the previous model
        B_q = theta(1:m,j);
        A_q = [1,theta(m+1:m+n, j)'];
        G_qu = (sum(B_q))/(sum(A_q));

        %%Dissimilarity measure for the gain of the models
        ar_pq(i,j) = abs(theta(end,i)-theta(end,j));
        Kappa_pq(i,j) = abs(G_pu - G_qu); %Dissimilarity measure

    end
end

%Check conditions
merge_gain = (Kappa_pq < kappa_K);
merge_affine = (ar_pq < kappa_r);
merge_olap = olap < kappa_olap;
merge_check = olap.*(merge_gain.*merge_affine.*merge_olap);
merge_check(merge_check <= 0) = NaN;

%Find minimum overlapping clusters
[min_col_vals, row_indices] = min(merge_check);
[min_olap, col] = min(min_col_vals);
row = row_indices(col);
im = row; jm = col;

%Dispaly clusters during the merging mechanism with the
% overlapping candidates
if enable_display_clusters
    subplot(32,1,1:19)
    plot([mi(im,1), mi(jm,1)],[mi(im,2),mi(jm,2)],'Color','r','linewidth',1.2);
end

%Merge
if any(merge_check(:)>0) %Check dissimilarity margin

    %Debugging display
    if merge_olap
        disp('Cluster Merge overlap')
    else
        disp('Cluster Merge transfer function')
    end

    %Merge clusters p and q
    n_p = floor(n_c(im));
    n_q = floor(n_c(jm));
    mi_p = mi(im,:);
    mi_q = mi(jm,:);
    n_pq = (n_p + n_q);
    mi(jm,:) = (n_p*mi_p + n_q*mi_q)/n_pq;
    n_c(jm,1) = n_pq;

    %Merge S_c
    E_p = ones(n_p,n_z);
    M_p = diag(mi_p);
    ZpTZp = (n_p-1)*(1/n_p)*S_c(:,:,im) + M_p'*(E_p')*E_p*M_p;

    E_q = ones(n_q,n_z);
    M_q = diag(mi_q);
    ZqTZq = (n_q-1)*(1/n_q)*S_c(:,:,jm) + M_q'*(E_q')*E_q*M_q;

    E_pq = ones(n_pq,n_z);
    M_pq = diag(mi(jm,:));
    Sigma = (1/(n_pq-1))*(ZpTZp + ZqTZq - M_pq'*(E_pq')*E_pq*M_pq);

    S_c(:,:,jm) = n_pq*Sigma;

    %Merge local parameters
    theta_p = theta(:, im);
    theta_q = theta(:, jm);
    theta(:,jm) = (theta_p + theta_q)/2;

    MSE(jm,1) = (n_q*MSE(jm) + n_p*MSE(im))/n_pq;
    b_fuzzy(jm)  = (n_q*b_fuzzy(jm) + n_p*b_fuzzy(im))/n_pq; %(n_q*b_fuzzy(jm) + n_p*b_fuzzy(im))/n_pq;
    P(:,:,jm) = (n_q*P(:,:,jm) + n_p*P(:,:,im))/n_pq;


    %Compute projection
%     Sigma(:,:,jm) = S_c(:,:,jm)/n_c(jm);
%     [~,S,V] = svd(Sigma(:,:,jm));
%     S = diag(sqrt(S));
%     delta(:,:,jm) = max(abs(V(1,:)'.*S));
    phi_c(jm,:) = (n_q*phi_c(jm,:) + n_p*phi_c(jm,:))/n_pq;

    %Remove one of the clusters
    S_c(:,:,im) = [];
    mi(im,:) = [];
    n_c(im) = [];
    b_fuzzy(im) = [];
    MSE(im) = [];
    theta(:,im) = [];
    P(:,:,im) = [];

%     Sigma(:,:,im) = [];
%     delta(:,:,im) = [];
    phi_c(im,:) = [];

    c = c - 1; %Remove one from the number of clusters

    merging = 1;
else
    merging = 0;
end