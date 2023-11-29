%Simulation of neuro-fuzzy models in parallel

%Load data
%valid = load('measurement_steps_all_2022_06_11_17_30.mat');
valid = load('measurement_steps_all_2022_06_14_11_28.mat');
%valid = load('measurement_steps_all_2022_06_16_21_40.mat');

ts = 2;
u_v = (valid.u -valid.u_min)/(valid.u_max -valid.u_min)*(u_max-u_min) + u_min;
x_v = valid.x;
k_0 = n+1;
x_v = x_v(find(diff(u_v),1)-n+1:end,:);
u_v = u_v(find(diff(u_v),1)-n+1:end);
y_v = x_v(:,2);
d_v = x_v(:,1);

%Add a fault of the output sensor
if enable_synthetic_fault
    y_v(6001:7000) = y_v(6001:7000) - 10;
end

z_v = [u_v, y_v, d_v];

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

%Signal initialization
u_est = u_v;
y_est = zeros(length(u_est),1);
y_est(1:(k_0+k_ahead-1)) = z_v(1:(k_0 + k_ahead-1),2);
y_dev = zeros(length(u_est),1);
y_est_top = zeros(length(u_est),1);
y_est_bot = zeros(length(u_est),1);
e_sim = zeros(length(u_est),1);
check_fault = zeros(length(u_est),1);
check_fault_f = zeros(length(u_est),1);
derivative_signal= zeros(length(u_est),1);
check_fault_f_der= zeros(length(u_est),1);
phi_P_phi_pred = zeros(length(u_est),1);
fault_diagnosis = zeros(length(u_est),n_phi);
phi_dist = zeros(length(u_est),n_phi);

y_pred(1:k_0) = z_v(1:k_0,2);
for k_v = (k_ahead + k_0):1:length(u_est)

    %Prediction setup
    u_pred = z_v((k_v-k_ahead-k_0+1):1:k_v,1);
    y_pred(1:k_0) = z_v((k_v-k_ahead-k_0+1):1:(k_v-k_ahead),2);
    d_pred = z_v((k_v-k_ahead-k_0+1):1:k_v,3);

    if 1

        %For all clusters
        for i = 1:1:c

            %Compute mmembership functions
            Gamma(i) = exp(-((z_v(k_v-1,1:2) - mi(i,1:2))*...
                pinv(delta_sim(1:2, 1:2, i))*(z_v(k_v-1,1:2) - mi(i,1:2))')^etta);

            %                Gamma(i) = exp(-((z_v(k_v-1,2) - mi(i,2))*...
            %                   pinv(delta_sim(2, 2, i))*(z_v(k_v-1,2) - mi(i,2))'));


            %Agregate the parameter vectors
            Theta(:,i) = theta(:,i)*Gamma(i);
            P_sum(:,:,i) = P(:,:,i)*Gamma(i);
            sigma(:,i) = MSE(i)*Gamma(i);
            Phi_c(i,:) = phi_c(i,:)*Gamma(i);

        end

        Theta = sum(Theta/sum(Gamma),2); %Normalize the parameter vector
        P_sum = sum(P_sum/sum(Gamma),3);
        sigma = sum(sigma/sum(Gamma));
        Phi_c = sum(Phi_c/sum(Gamma),1);
    end

    %Compute k steps ahead prediction
    for ks = (k_0):1:(k_ahead+k_0)

        %
        if 0
            %For all clusters
            for i = 1:1:c

                %Compute mmembership functions
                Gamma(i) = exp(-(([u_pred(ks-1), y_pred(ks-1)] - mi(i,1:2))*...
                    pinv(delta_sim(1:2, 1:2, i))*([u_pred(ks-1), y_pred(ks-1)] - mi(i,1:2))')^etta);

                %Agregate the parameter vectors
                Theta(:,i) = theta(:,i)*Gamma(i);

                if ks == (k_ahead+k_0) %Only the last confidence interval is used
                    P_sum(:,:,i) = P(:,:,i)*Gamma(i);
                    sigma(:,i) = MSE(i)*Gamma(i);
                    Phi_c(i,:) = phi_c(i,:)*Gamma(i);
                end

            end

            Theta = sum(Theta/sum(Gamma),2); %Normalize the parameter vector
        end

        %Regressor for the consequence part
        Phi = [u_pred(ks-1:-1:ks-m)', -y_pred(ks-1:-1:ks-n), d_pred(ks-1:-1:ks-m)', 1]';
        y_pred(ks) = Phi'*Theta;

    end

    %Save the last output of the prediction and prediction intervals
    y_est(k_v) = y_pred(end);
    y_dev(k_v,1) = sqrt(sigma*(1 + Phi'*P_sum*Phi));
    phi_P_phi_pred(k_v) = Phi'*P_sum*Phi;
    y_est_top(k_v,1) = y_est(k_v) + t_deviations*y_dev(k_v);
    y_est_bot(k_v,1) = y_est(k_v) - t_deviations*y_dev(k_v);
    e_sim(k_v) = mean((z_v((k_v-k_ahead-k_0+1):1:k_v,2)-y_pred').^2);

    fault_diagnosis(k_v,:) = abs((((P_sum*Phi)/(Phi'*P_sum*Phi+1))./Theta)');
    phi_dist(k_v,:) =  abs((Phi_c-Phi')'.*P_sum*(Phi_c-Phi')'./Theta);

    A_f = [1, Theta(m+1:m+n)]';
    check_fault(k_v)  = abs(max(phi_dist(k_v,:))*(e_sim(k_v)*(phi_P_phi_pred(k_v))));
    check_fault_f(k_v) = sum(A_f)*check_fault(k_v) - A_f(2:(n+1))*check_fault_f(k_v-1:-1:k_v-n);

end

%Model error
e_sim_sum = sum(e_sim.^2)/length(y_est);
disp(['N-steps-ahead prediction MSE = ', num2str(e_sim_sum)])

%Display prediction
figure(3); subplot(30,1,1:18);
set(gcf, 'Position', 1.1*[100 100 800 600])
hold off;
k_0_plot = 4;
h_y_v = plot(y_est(k_0_plot:end),'r'); hold on;
h_z_v_deviation = plot(y_est_top(k_0_plot:end),'k-.');
plot(y_est_bot(k_0_plot:end),'k-.');
h_z_v_y = plot(z_v(k_0_plot:end,2),'b');
h_z_v_d = plot(z_v(k_0_plot:end,3),'m');
xlim('tight')
ylim([20,80])

%Display true fault
k_steps = 0:1:length(u_est);
if enable_synthetic_fault
    h_fault_start = area(80*( ((k_steps>2150).*(k_steps<4900)) + ...
        ((k_steps>8134).*(k_steps<9640)) + ...
        ((k_steps>6000).*(k_steps<7000))), ...
        'EdgeColor','none','FaceColor','g','FaceAlpha',0.1); hold on;
else
    h_fault_start = area(80*( ((k_steps>2150).*(k_steps<4900)) + ...
        ((k_steps>8134).*(k_steps<9640))), ...
        'EdgeColor','none','FaceColor','g','FaceAlpha',0.1); hold on;
end

%Display detected fault
h_fault = area(80*(check_fault_f > fault_threshold), ...
    'EdgeColor','none','FaceColor','r','FaceAlpha',0.2);

%Figure settings
ax = gca;
set(gca,'fontname','Helvetica', 'FontSize', 12);
h=legend([h_y_v,h_z_v_deviation,h_z_v_y,h_z_v_d,h_fault,h_fault_start],...
    'Model output $\hat{y}$',...
    'Prediction interval $\hat{y}  \pm t \sqrt{cov(y-\hat{y})}$',...
    'System output measurement  $y$',...
    'Reservoir hot water $d$',...
    'Fault detected',...
    'Actual fault',...
    'location','southeast','Interpreter','latex');
set(h,'FontSize',10)
axis manual
xlim tight
xlabel('Time step $k$')
ylabel('Output signal')
%title('Fault detection')
set(h,'FontSize',10)
set(gca,'FontSize', 12);
name = ['heat_exchanger_prediction_intervals_prediction_interval_fault','.pdf'];
exportgraphics(gcf,name,'BackgroundColor','none');

figure(4)
subplot(2,1,1); hold off;
h_y_dev_top = plot(t_deviations*y_dev(k_0_plot:end) ,'k--'); hold on;
h_y_dev_bot = plot(-t_deviations*y_dev(k_0_plot:end),'k--');
h_error = plot(e_sim,'b');
ylabel('Prediction error')
ylim("tight");xlim("tight");

h = legend([h_y_dev_top,h_error],'Prediction intervals','Prediction error','Interpreter','latex');

subplot(2,1,2); hold off;
if enable_synthetic_fault
    h_fault_start = area(max(max(fault_diagnosis))*( ((k_steps>2150).*(k_steps<4900)) + ...
        ((k_steps>8134).*(k_steps<9640)) + ...
        ((k_steps>6000).*(k_steps<7000))), ...
        'EdgeColor','none','FaceColor','g','FaceAlpha',0.1); hold on;
else
    h_fault_start = area(max(max(fault_diagnosis))*( ((k_steps>2150).*(k_steps<4900)) + ...
        ((k_steps>8134).*(k_steps<9640))), ...
        'EdgeColor','none','FaceColor','g','FaceAlpha',0.1); hold on;
end
h_fault = area(max(check_fault_f)*(check_fault_f > fault_threshold), ...
    'EdgeColor','none','FaceColor','r','FaceAlpha',0.3);
set(gca,'fontname','Helvetica', 'FontSize', 12); hold on;
plot(check_fault,'b');
plot(check_fault_f,'m--');
legend('Fault','Detected fault','Alarm signal','Filtered alarm signal','location','best')
xlim("tight");ylim([0,max(check_fault_f)]);
xlabel('Time step $k$')
ylabel('Alarm signal')

figure(5)
subplot(2,1,1); hold off;
if enable_synthetic_fault
    h_fault_start = area(max(max(fault_diagnosis))*( ((k_steps>2150).*(k_steps<4900)) + ...
        ((k_steps>8134).*(k_steps<9640)) + ...
        ((k_steps>6000).*(k_steps<7000))), ...
        'EdgeColor','none','FaceColor','g','FaceAlpha',0.1); hold on;
else
    h_fault_start = area(max(max(fault_diagnosis))*( ((k_steps>2150).*(k_steps<4900)) + ...
        ((k_steps>8134).*(k_steps<9640))), ...
        'EdgeColor','none','FaceColor','g','FaceAlpha',0.1); hold on;
end
plot(fault_diagnosis);
legend('Fault','Input signal', 'Output signal', 'Reservoir signal', 'Affine constant')
ylim("tight");xlim("tight");
ylabel('Kalman gain')

subplot(2,1,2); hold off;
if enable_synthetic_fault
    h_fault_start = area(max(max(phi_dist))*( ((k_steps>2150).*(k_steps<4900)) + ...
        ((k_steps>8134).*(k_steps<9640)) + ...
        ((k_steps>6000).*(k_steps<7000))), ...
        'EdgeColor','none','FaceColor','g','FaceAlpha',0.1); hold on;
else
    h_fault_start = area(max(max(phi_dist))*( ((k_steps>2150).*(k_steps<4900)) + ...
        ((k_steps>8134).*(k_steps<9640))), ...
        'EdgeColor','none','FaceColor','g','FaceAlpha',0.1); hold on;
end
plot(phi_dist);
legend('Fault','Input signal', 'Output signal', 'Reservoir signal', 'Affine constant')
ylim("tight");xlim("tight");
ylabel('Distance to regressor')
xlabel('Time step $k$')