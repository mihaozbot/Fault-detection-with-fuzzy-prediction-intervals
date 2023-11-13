valid = load('measurement_steps_all_2022_06_11_17_30.mat');
valid = load('measurement_steps_all_2022_06_14_11_28.mat');
%valid = load('measurement_steps_all_2022_06_16_21_40.mat');
t_deviations = 1;
ts = 2;
u_v = (valid.u -valid.u_min)/(valid.u_max -valid.u_min)*(u_max-u_min) + u_min;
x_v = valid.x;
k_0 = n+1;
x_v = x_v(find(diff(u_v),1)-n+1:end,:);
u_v = u_v(find(diff(u_v),1)-n+1:end);

%Simulation of neuro-fuzzy models in parallel
paren = @(x, varargin) x(varargin{:});

y_v = x_v(:,2);
d_v = x_v(:,1);
z_v = [u_v,y_v, d_v];

u_est = u_v;
y_est = y_v(1:k_0);

%Calculation of clusters projection parameter
delta_sim = [];
for i = 1:1:c
    Sigma(:,:,i) = S_c(:,:,i)/n_c(i);
    [~,S,V] = svd(Sigma(:,:,i));
    S = diag(sqrt(S));
    delta_sim(:,:,i) = diag([max(abs(V(1,:)'.*S)),max(abs(V(2,:)'.*S))]);
end

%Simulation
Gamma = [];
sigma = [];
Theta = [];

I_e_f = zeros(k_0,1);

z_f = z_v(1:(k_0-1),2);
e_sim = zeros(1:(k_0-1),1);
e_f  = e_sim;

for k_v = k_0:1:length(u_est)

    for i = 1:1:c
        Gamma(i) = exp(-(([u_est(k_v-1),y_est(k_v-1)] - mi(i,1:2))*...
            inv(delta_sim(1:2,1:2,i))*([u_est(k_v-1),y_est(k_v-1)] - mi(i,1:2))')^etta);
        P_sum(:,:,i) = P(:,:,i)*Gamma(i);
        Theta(:,i) = theta(:,i)*Gamma(i);
        sigma(:,i) = MSE(i)*Gamma(i);

    end

    P_sum = sum(P_sum/sum(Gamma),3);
    Theta = sum(Theta/sum(Gamma),2);
    sigma = sum(sigma/sum(Gamma));

    for i = 1:c

    end

    z_f(k_v) = sum(A_f)*z_v(k_v,2) - A_f(2:(n+1))*z_f(k_v-1:-1:k_v-n);
    
    %Regressor
    if enable_disturbance
        Phi = [u_est(k_v-1:-1:k_v-m)', -y_est(k_v-1:-1:k_v-n), z_v(1,3), 1]';
    elseif ~enable_disturbance
        Phi = [u_est(k_v-p-1:-1:k_v-m)', -y_est(k_v-1:-1:k_v-n), 1]';
    end
    
    A_f = [1,Theta(m+1:m+n)']';

    y_est(k_v) = Phi'*Theta;
    y_dev(k_v) = sqrt(sigma*(1+Phi'*P_sum*Phi));
    y_est_top(k_v) = Phi'*Theta  + t_deviations*y_dev(k_v);
    y_est_bot(k_v) = Phi'*Theta  - t_deviations*y_dev(k_v);
    e_sim(k_v) = (z_v(k_v,2)-y_est(k_v));
    e_f(k_v) = sum(A_f)*e_sim(k_v) - A_f(2:(n+1))*e_f(k_v-1:-1:k_v-n);
    d_e_f(k_v) = (e_f(k_v)-e_f(k_v-1))/ts;
    I_e_f(k_v) = I_e_f(k_v-1) + ts*e_f(k_v);
 
end

%Model error
e_sim_sum = sum(e_sim.^2)/length(y_est);
disp(['e_sim_sum = ', num2str(e_sim_sum)])

figure(3); subplot(30,1,1:18);
hold off;
k_0_plot = 4;
h_y_v = plot(y_est(k_0_plot:end),'r'); hold on;
h_z_v_deviation = plot(y_est_top(3:end),'k-.');
plot(y_est_bot(k_0_plot:end),'k-.');
h_z_v_y = plot(z_v(k_0_plot:end,2),'b'); 
h_z_v_d = plot(z_v(k_0_plot:end,3),'m');  

xlim('tight')
ylim([20,80])

if synthetic_fault
    h_fault_start = area(max(max(phi_dist))*( ((k_steps>2150).*(k_steps<4900)) + ...
        ((k_steps>8134).*(k_steps<9640)) + ...
        ((k_steps>6000).*(k_steps<7000))), ...
        'EdgeColor','none','FaceColor','g','FaceAlpha',0.1); hold on;
else
    h_fault_start = area(max(max(phi_dist))*( ((k_steps>2150).*(k_steps<4900)) + ...
        ((k_steps>8134).*(k_steps<9640))), ...
        'EdgeColor','none','FaceColor','g','FaceAlpha',0.1); hold on;
end

axis manual
%plot(sigma_les_est,'r--')
xlabel('Step [$k$]')
ylabel('Output signal')

ax = gca;
ax.Toolbar.Visible = 'off';
%set(ax,'fontname','Times', 'FontSize', 14);
%set(ax,'fontname','cmr12', 'FontSize', 14);
set(gca,'fontname','Helvetica', 'FontSize', 14);
h=legend([h_y_v,h_z_v_deviation,h_z_v_y,h_z_v_d,h_fault,h_fault_start],...
    'Model output $\hat{y}$',...
    'Prediction interval $\hat{y}  \pm t \sqrt{cov(y-\hat{y})}$',...
    'System output measurement  $y$',...
    'Reservoir hot water $d$',...
    'Fault detected',...
    'Actual fault',...
    'location','southeast','Interpreter','latex');
set(h,'FontSize',14)
h = title('Fault detection');
%set(h,'FontSize',14,'FontWeight','Normal','fontname','Helvetica')

subplot(30,1,20:30); hold off;
h_e= plot(e_sim(k_0_plot:end),'r'); hold on;
h_e_f = plot(e_f(k_0_plot:end),'b'); hold on;
h_y_dev_top = plot(t_deviations*y_dev(k_0_plot:end),'k--'); 
h_y_dev_bot = plot(-t_deviations*y_dev(k_0_plot:end),'k--'); 
set(gca,'fontname','Helvetica', 'FontSize', 14);

name = ['prediction_interval_fault','.pdf'];
saveas(gcf, name);

figure(4)
subplot(16,1,1:5); hold off;
plot(e_sim,'r'); hold on;
plot(e_f,'b')
ylabel('Filtered error $e_f$')
ylim('tight')
xlim('tight')
title('Filtered simulation error')
set(gca,'fontname','Helvetica', 'FontSize', 14);
legend('Prediction error $e(k)$','Filtered error $e_f(k) = \frac{A_f(1)}{A_f(q^{-1})}e(k)$','location','best','Interpreter','latex');

%set(gca,'xticklabel',[])
subplot(16,1,7:10)
plot(I_e_f,'b')
ylabel('Integral of $e_f$')
ylim('tight')
xlim('tight')
set(gca,'fontname','Helvetica', 'FontSize', 14);
%set(gca,'xticklabel',[])

subplot(16,1,12:16)
plot(d_e_f,'b')
xlabel('Step k')
ylabel('Derivative of $e_f$')
ylim('tight')
xlim('tight')
set(gca,'fontname','Helvetica', 'FontSize', 14);
