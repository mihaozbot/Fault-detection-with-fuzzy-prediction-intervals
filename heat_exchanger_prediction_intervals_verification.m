%Simulation of neuro-fuzzy models in parallel

k_0 = n+1; %First time step for the simulation based on the order of the sys

%Load data
x_ver = x; 
u_ver = u;
paren = @(x, varargin) x(varargin{:});

%Setup data
y_ver = x_ver(:,2);
d_ver = x_ver(:,1);
z_ver = [u_ver, y_ver, d_ver];

%Initialize model prediction
u_est_ver = u_ver;
y_est_ver = zeros(size(u_est_ver));
y_est_ver(1:k_0) = y_ver(1:k_0);
y_est_ver_top = zeros(size(u_est_ver));
y_est_ver_bot = zeros(size(u_est_ver));

%Calculation of clusters projection parameter
Sigma = zeros(size(S_c));
delta_sim = zeros(2, 2, c);
for i = 1:1:c
    Sigma(:,:,i) = S_c(:,:,i)/n_c(i);
    [~,S,V] = svd(Sigma(:,:,i));
    S = diag(sqrt(S));
    delta_sim(:,:,i) = diag([max(abs(V(1,:)'.*S)),max(abs(V(2,:)'.*S))]);
end

%Simulation
for k_v = k_0:1:length(u_est_ver)

    %Compute the activation
    Gamma = zeros(c, 1);
    P_sum_act = zeros(size(P));
    Theta_act = zeros(size(theta));
    sigma_act = zeros(size(MSE));

    for i = 1:1:c
         Gamma(i) = exp(-(([u_est_ver(k_v-1),y_est_ver(k_v-1)] - mi(i,1:2))*...
             inv(delta_sim(1:2,1:2,i))*([u_est_ver(k_v-1),y_est_ver(k_v-1)] - mi(i,1:2))')^etta);
%                 Gamma(i) = exp(-(([u_est_ver(k_v-1)] - mi(i,1))*...
%             inv(delta_sim(1,1,i))*([u_est_ver(k_v-1)] - mi(i,1))')^etta);
        P_sum_act(:,:,i) = P(:,:,i)*Gamma(i);
        Theta_act(:,i) = theta(:,i)*Gamma(i);
        sigma_act(i,1) = MSE(i).*Gamma(i);
    end

    %Normalize the parameters based on the activation
    P_sum = sum(P_sum_act/sum(Gamma),3);
    Theta = sum(Theta_act/sum(Gamma),2);
    sigma = sum(sigma_act/sum(Gamma));

    %Regressor
    if enable_disturbance %Use the hot water reservoar in the regressor
        Phi = [u_est_ver(k_v-1:-1:k_v-m)', -y_est_ver(k_v-1:-1:k_v-n), z_ver(k_v,3), 1]';
    else
        Phi = [u_est_ver(k_v-1:-1:k_v-m)', -y_est_ver(k_v-1:-1:k_v-n), 1]';
    end

    %Output and confidence interval
    y_est_ver(k_v) = Phi'*Theta;
    y_dev_ver = sqrt(sigma*(1+Phi'*P_sum*Phi));
    y_est_ver_top(k_v) = Phi'*Theta  + t_deviations*y_dev_ver;
    y_est_ver_bot(k_v) = Phi'*Theta  - t_deviations*y_dev_ver;

end

%Model error
e_sim  = (z_ver(:,2)-y_est_ver).^2;
e_sim_sum = sum(e_sim)/length(y_est_ver);
disp(['Simulation MSE: ', num2str(e_sim_sum)])

%Plot results of the verification
figure(2); hold off;
subplot(2,1,1); hold off;
plot(z_ver(k_0:end,2),'b'); hold on;
plot(y_est_ver(k_0:end),'r')

plot(y_est_ver_top(k_0:end),'k')
plot(y_est_ver_bot(k_0:end),'k')

xlabel('Time step $k$')
ylabel('Output signal [°C]')
xlim('tight'),ylim('tight')
ax = gca;
ax.Toolbar.Visible = 'off';
set(ax,'fontname','Times', 'FontSize', 10);
h=legend('System output $y$','Model output $\hat{y}$','Confidence interval',...
    'location','best','Interpreter','latex');
set(h,'FontSize',10)
%h = title('Simulation in parallel');
set(h,'FontSize',10,'FontWeight','Normal','fontname','Times')
