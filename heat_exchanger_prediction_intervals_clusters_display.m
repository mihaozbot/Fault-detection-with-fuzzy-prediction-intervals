set(groot,'defaultLineLineWidth',1.5)

%Select distinct colors for the clusters
if c > 1
    color = pmkmp(c,'IsoL');
else
    color = [1,0,0];
end

%Calculation of clusters projection parameter
delta_plot = [];
Sigma_plot = [];
for i = 1:1:c
    Sigma_plot(:,:,i) = S_c(:,:,i)/n_c(i);
    [~,S,V] = svd(Sigma_plot(:,:,i));
    S = diag(sqrt(S));
    Vu = V(1,:)';
    Vy = V(2,:)';
    delta_plot(:,:,i) = diag([max(abs(Vu.*S)),max(abs(Vy.*S))]);
    n_pts = 100; % Number of points around ellipse
    kot = 0:pi/n_pts:2*pi; % angles around a circle
    [eigvec, eigval] = eig(Sigma_plot([1:2],[1:2],i));
    xy = [cos(kot'),sin(kot')] *sqrt(eigval) * eigvec';
    x_sigma(:,i) = xy(:,1);
    y_sigma (:,i)= xy(:,2);
end

%Compute the normalized membership functions
u_cluster = u_min:0.01:u_max;
Gamma_norm = zeros(length(u_cluster),c);
for k_disp = 1:1:length(u_cluster)

    %Membership functions
    Gamma = zeros(c,1);
    for i = 1:1:c
        Gamma(i) = exp(-(1/delta_plot(1,1,i))*((u_cluster(k_disp) - mi(i,1))*(u_cluster(k_disp) - mi(i,1)))^etta);
    end
    Gamma_norm(k_disp,:) = Gamma./sum(Gamma);
end

%Plot clusters 
figure(1);
subplot(32,1,1:19); hold off;
plot(u(1:k),y(1:k),'k.'); hold on;
axis manual
for i = 1:1:c
    plot(mi(i,1),mi(i,2),'o','Color',color(i,:),'linewidth',1.2)
    plot(x_sigma(:,i) + mi(i,1),y_sigma(:,i)+ mi(i,2),'Color',color(i,:))
    xy = [cos(kot'),sin(kot')] *delta_plot(:,:,i);
    plot(mi(i,1)+xy(:,1),mi(i,2)+xy(:,2),'--','Color',color(i,:));
end
xlim([u_min,u_max])
xlabel('Input signal $u$')
ylabel('System output $y$')
ax = gca;
ax.Toolbar.Visible = 'off';
set(ax,'fontname','Times', 'FontSize', 14);

K_c = NaN(c,m+m);
N_c = NaN(c,1);
for i = 1:1:c

        K_c(i,:) = [(theta(1,i)),(theta(3,i))]/(sum([1,theta((m+1):(m+n), i)']));
        N_c(i) = theta(end, i)/(sum([1,theta((m+1):(m+n),i)']));

        d_u = 0.1;
        p3 = plot(u_min:(d_u):u_max,K_c(i,1)*((u_min:(d_u):u_max)) + K_c(i,2)*mi(c,3) + N_c(i),'--','Color',color(i,:),'linewidth',1);

end

%Plot normalized activation functions
subplot(32,1,25:32); hold off;
for i = 1:1:c
    plot(u_cluster,Gamma_norm(:,i),'Color',color(i,:),'linewidth',1.2); hold on;
end
xlabel('Input signal $u$')
ylabel('Cluster memb. $\rho_i$')
yticks([0,1])
ax = gca;
ax.Toolbar.Visible = 'off';
set(ax,'fontname','Times', 'FontSize', 14);


