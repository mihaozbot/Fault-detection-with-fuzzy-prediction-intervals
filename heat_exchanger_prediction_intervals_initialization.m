
%Load data
data = load('measurement_steps_all_2022_06_11_17_30.mat');

k_0 = n+1+1; %First time step of the identificaiton based on the order of the system

%Scale data to correct range
u = (data.u -data.u_min)/(data.u_max -data.u_min)*(u_max-u_min) + u_min;
x = data.x;
x = x(find(diff(u),1)*2-n+1:end,:);
u = u(find(diff(u),1)*2-n+1:end);

%Reshape data for the model
u = u(:);
y = x(:,2);
d = x(:,1);
z = [u(1:k_0), y(1:k_0), d(1:k_0)];
d_f = z(:,3);

%Initialize the filtered signals
u_f = z(:,1);
y_f = z(:,2);

%Setup the sytem dimensions
if enable_disturbance
    n_phi = m+n+1+1;
    n_z = 3;
else
    n_phi = m+n+1;
    n_z = 2;
end

%Initialize the information matrix
P0 = P_alpha*eye(n_phi);

%Initialize the confidence interval variables
MSE = 0; %Mean squared error
b_fuzzy = 0; %Fuzzy membership acumulated
phi_P_phi = 1*ones(k_0,1); %Confidence interval (core or without the variance)
%y_m_w = z(:,2);
y_m_wo = z(:,2); %Model estimation without the last added rule
%e_w = zeros(k_0,1);
e_wo = zeros(k_0,1); %Model error

%Add and initialize new cluster
c = 0; %Increase cluster counter

%This is the whole model representation. The first rules is added inside the
%identificaiton loop
% P(:,:,c) = P0;
% n_c(c) = 1; %Set number of samples of cluster to 1
% mi(c,:) = z(end-1,:);
% S_c(:,:,c) = zeros(n_z,n_z);
% theta(:,c) = zeros(n_phi,1);
