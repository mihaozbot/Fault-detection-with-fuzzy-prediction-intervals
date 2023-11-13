

%Add and initialize new rule
c = c + 1; %Increase cluster counter

%Reset the model prediction
%y_m_w(k:-1:k-n) = z(k:-1:k-n,2);
y_m_wo(k:-1:k-n) = z(k:-1:k-n,2);
y_m(k:-1:k-n) = z(k:-1:k-n,2);

%Initialize new rule
n_c(c,1) = 1; %Set number of samples of cluster to 1
mi(c,:) = z(end,1:n_z); %Cluster mean
S_c(:,:,c) = S_0*eye(n_z,n_z); %Clusters covariance matrix
theta(:,c) = zeros(n_phi,1); %Consequence parameters
P(:,:,c) = P0; %Inforamtion or covariance matrix for param. optimization
MSE(c,1) = 0; %Mean squared error
b_fuzzy(c,1) = 1; %Accumulated activation
phi_P_phi(k) = 1; %Confidence interval
Gamma(c,1) = 1; %Membership activation
jg = c;
phi_c(c,:) = phi;
