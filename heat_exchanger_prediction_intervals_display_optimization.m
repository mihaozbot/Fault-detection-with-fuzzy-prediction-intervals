set(0,'defaultTextInterpreter','latex')
set(groot,'defaultLineLineWidth',1.5)

paren = @(x, varargin) x(varargin{:});


%Calculation of the covariance matrices
if n_z ==2
for i = 1:1:c

    Sigma(:,:,i) = (S_c(:,:,i)./n_c(i));
    n_pts = 100; % Number of points around ellipse
    kot = 0:pi/n_pts:2*pi; % angles around a circle
    [eigvec, eigval] = eig(Sigma(:,:,i)); % Compute eigen-stuff
    xy = [cos(kot'),sin(kot')] *sqrt(eigval) * eigvec'; % Transformation
    x_sigma(:,i) = xy(:,1);
    y_sigma(:,i) = xy(:,2);

end
end
%Plot color vector

%color = distinguishable_colors(c);
%LinearL,
color = pmkmp(c,'IsoL');
% for i = 1:1:length(p_all_unsorted)
%     [~,i_closest(i)] = min(abs(mi(:,1)-p_all_unsorted(i)));
% end
% 

% 
% figure(3)
% subplot(2,1,1);hold off
% %plot(w(1:k),'b'); hold on;
% plot(z(:,2),'r'); axis manual;
% plot(y_m_w,'k')
% plot(y_m_wo,'k--')
% ylabel('Output signals')
% legend('Reference signal $r$','System output $y$','Current model output $\hat{y}_c$','Whole model output $\hat{y}$','Location','best','Interpreter','Latex')
% ax = gca;
% ax.Toolbar.Visible = 'off';
% set(ax,'fontname','Times', 'FontSize', 14);
% subplot(4,1,3);
% plot(u,'b');hold on
% ylabel(["Control";"signal $u$"])
% yticks(linspace(min(u),max(u),3))
% ytickformat('%,.0f')
% ax = gca;
% ax.Toolbar.Visible = 'off';
% set(ax,'fontname','Times', 'FontSize', 14);
% % subplot(4,1,4)
% % %plot(e_p,'r')
% % ylabel(["Control";" error $e_c$"])
% % xlabel('Time step')
% % yticks(linspace(min(e_p),max(e_p),3))
% % ytickformat('%,.1f')
% %legend('u','e','gamma_c','Location','best')
% pause(0.01)
% ax = gca;
% ax.Toolbar.Visible = 'off';
% set(ax,'fontname','Times', 'FontSize', 14);

% figure
% plot(Theta_all);hold on;
% plot(phi_P_phi_posteriori > kappa_c)
% plot(z(:,2),'r')

figure;
%plot(50*(stability_and_confidence)); hold on;
plot(z(:,2),'k')

% 
% figure
% subplot(2,1,1);hold off
% plot(w(1:k),'b'); hold on;
% plot(z(:,2),'r'); axis manual;


figure;
subplot(2,1,1)
plot(y); hold on;
ylabel('Signal [ °C]')
xlabel('Time step $k$')
ax = gca;
ax.Toolbar.Visible = 'off';
set(ax,'fontname','Times', 'FontSize', 14);
h = legend('System output $T_{sp}$','Reservoir signal $T_{ep}$','location','northeast','Interpreter','latex');
set(h,'FontSize',14)

subplot(2,1,2)
plot(u,'b')
h = title('Potek vzbujalnega signala');
set(h,'FontSize',14,'FontWeight','Normal','fontname','Times')
ylabel('Input signal $u$')
xlabel('Time step $k$')
ylim([2,22])
ax = gca;
ax.Toolbar.Visible = 'off';
set(ax,'fontname','Times', 'FontSize', 14);

figure;
subplot(8,1,1:5)
 hold on;

plot(y,'b-','linewidth',1.2);
plot(y_f,'r--','linewidth',1.2)
%plot(r_all);
ylabel('Output signal')
xlabel('Time step $k$')
ax = gca;
ax.Toolbar.Visible = 'off';
set(ax,'fontname','Times', 'FontSize', 14);
h = legend('System output $y$','Filtered system output $y_f$','location','northeast','Interpreter','latex');
set(h,'FontSize',14)
%h = title('Filtriran izhodni signal');
set(h,'FontSize',14,'FontWeight','Normal','fontname','Times')
%ylim([-3,4.2])
subplot(8,1,7:8)
plot(u,'b','linewidth',1.2)
%h = title('Potek vzbujalnega signala');
set(h,'FontSize',14,'FontWeight','Normal','fontname','Times')
ylabel('Input signal $u$')
xlabel('Time step $k$')
%ylim([2,22])
ax = gca;
ax.Toolbar.Visible = 'off';
set(ax,'fontname','Times', 'FontSize', 14);
name = ['FRLS_heuristika_filtered_signal','.pdf'];
exportgraphics(gcf, name,'ContentType','vector');

%Parameters estimation
if n == 2
    figure

    subplot(32,1,1:5);hold on
    plot(M.azz_1(u_all),'r','linewidth',1.2)
    plot(theta_check(2,:),'b','linewidth',1.2);hold on
    %plot(theta_check(2,:)+theta_sigma_all(2,:),'b--')
    %plot(theta_check(2,:)-theta_sigma_all(2,:),'b--')
    xlim([0,length(u_all)])
    ylabel('${a_1}$')
    xlabel('Time step $k$')
    %h = legend('Theoretical par.','Identified par.','location','northwest');
    ax = gca;
    ax.Toolbar.Visible = 'off';
    set(ax,'fontname','Times', 'FontSize', 14);

    %set(h,'FontSize',14,'FontWeight','Normal','fontname','Times')
    %set(h,'FontSize',14)

    subplot(32,1,11:15);hold on
    plot(M.azz_0(u_all),'r','linewidth',1.2)
    plot(theta_check(3,:),'b','linewidth',1.2)
    %plot(theta_check(3,:)+theta_sigma_all(3,:),'b--')
    %plot(theta_check(3,:)-theta_sigma_all(3,:),'b--')
    xlim([0,length(u_all)])
    ylabel('$a_2$')
    xlabel('Time step $k$')
    ax = gca;
    ax.Toolbar.Visible = 'off';
    set(ax,'fontname','Times', 'FontSize', 14);

    subplot('Position',[0.15,0.15,0.32,0.25]);hold on
    plot(M.bzz_1(u_all)+M.bzz_0(u_all),'r','linewidth',1.2)
    plot(theta_check(1,:),'b','linewidth',1.2)
    %plot(theta_check(1,:)+theta_sigma_all(1,:),'b--')
    %plot(theta_check(1,:)-theta_sigma_all(1,:),'b--')
    xlim([0,length(u_all)])
    ylabel('${b_1}$')
    ylim([0.005,0.02]);
    xlabel('Time step $k$')
    ax = gca;
    ax.Toolbar.Visible = 'off';
    set(ax,'fontname','Times', 'FontSize', 14);

    subplot('Position',[0.6,0.15,0.32,0.25]);hold on
    plot(M.rzz(u_all'),'r.','linewidth',1.2)
    plot(theta_check(4,:),'b','linewidth',1.2)
    %plot(theta_check(4,:)+theta_sigma_all(4,:),'b--')
    %plot(theta_check(4,:)-theta_sigma_all(4,:),'b--')
    ylabel('${r}$')
    ylim([-0.2,0.2]);
    xlabel('Time step $k$')
    xlim([0,length(u_all)])
    yticks([-0.2,-0.1,0,0.1,0.2])
    ax = gca;
    ax.Toolbar.Visible = 'off';
    set(ax,'fontname','Times', 'FontSize', 14);
    %h = sgtitle('Optimizacija parametrov modela');
    set(h,'FontSize',14,'FontWeight','Normal','fontname','Times')
    name = ['FRLS_heuristika_parametri','.pdf'];
    exportgraphics(gcf, name,'ContentType','vector');

elseif (n == 1)

    figure
    subplot(3,1,1);hold on
    %plot(M.bzz_1(u_all)+M.bzz_0(u_all),'r')
    plot(theta_check(1,:),'b')
    ylabel('${b_1}$')
    xlim([0,length(u)])
    ylim([-0.0,0.2]);
    xlabel('Time step $k$')
    ax = gca;
    ax.Toolbar.Visible = 'off';
    set(ax,'fontname','Times', 'FontSize', 14);

    subplot(3,1,2);hold on
    plot(theta_check(2,:),'b');hold on
    ylabel('${a_1}$')
    xlim([0,length(u)])
    ylim([-1,-0.9]);
    xlabel('Time step $k$')
    ax = gca;
    ax.Toolbar.Visible = 'off';
    set(ax,'fontname','Times', 'FontSize', 14);

    subplot(3,1,3);hold on
    plot(theta_check(3,:),'b')
    ylabel('$r$')
    xlim([0,length(u)])
    ylim([-1,2]);
    xlabel('Time step $k$')
    ax = gca;
    ax.Toolbar.Visible = 'off';
    set(ax,'fontname','Times', 'FontSize', 14);
    h = sgtitle('Optimizacija parametrov modela');
    set(h,'FontSize',14,'FontWeight','Normal','fontname','Times')
    %name = ['FRLS_PWL_drevo_optimizacija_parametrov_',num2str(r),'.pdf'];
end

if n_z ==2
    %Input-output space
    figure; hold on;
    title('Input-output space')
    plot(u_all,y_all,'b.')
    for i = 1:1:c
        h1 = plot(mi(i,1),mi(i,2),'o','markerfacecolor',color(i,:));
        plot(x_sigma(:,i) + mi(i,1),y_sigma(:,i)+ mi(i,2), 'Color',color(i,:))
    end
    %set(h1, 'markerfacecolor', get(h1, 'color'));
    xlabel('Process input $u$')
    ylabel('Process output $y$')
    
elseif n_z ==3

    figure; hold on;
    title('Input-output space')
    plot3(u,y,z(:,3),'b.')
    for i = 1:1:c
        h1 = plot_gaussian_ellipsoid(mi(i,:),Sigma(:,:,i),1,10);
        set(h1,'FaceColor',color(i,:), ...
            'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none')
        plot3(mi(i,1),mi(i,2),mi(i,3),'ok');
        %plot3(x_sigma(:,i) + mi(i,1),y_sigma(:,i)+ mi(i,2),mi(i,3)*ones(length(y_sigma),1), 'Color',color(i,:))
    end
    %set(h1, 'markerfacecolor', get(h1, 'color'));
    xlabel('Process input $u$')
    ylabel('Process output $y$')
    zlabel('Reference temperature $z$')
    %zlim([85,95])
    view(75,15)
end
% if evolve_check
% %Sistem gain
% figure; hold on;
% plot(G_1,'b')
% 
% %plot(1:1:length(G_1),1./(1+u.^2),'r')
% ylabel('Enosmerno ojačenje $K$')
% ylim([-0,6])
% xlabel('Time step $k$')
% ax = gca;
% ax.Toolbar.Visible = 'off';
% set(ax,'fontname','Times', 'FontSize', 14);
% h = legend('Teoretična vrednost ojačenja','Identificirana vrednost ojačenja','location','southeast');
% set(h,'FontSize',14,'FontWeight','Normal','fontname','Times')
% h = title('Enosmerno ojačenje');
% set(h,'FontSize',14,'FontWeight','Normal','fontname','Times')
% end
% %
% [~, i_last_unique] = unique(u_all,'last');
% [~, i_max_unique] = max(i_last_unique);
% i_last_unique(i_max_unique) = [];
% i_last_unique_sort = sort(i_last_unique);
% u_last_unique = u_all(i_last_unique_sort);
% y_last_unique = y_all(i_last_unique_sort);
% 
% [~, i_first_unique] = unique(u_all,'first');
% [~, i_min_unique] = min(i_first_unique);
% i_first_unique(i_min_unique) = [];
% i_first_unique_sort = sort(i_first_unique);
% u_first_unique = u_all(i_first_unique_sort);
% y_first_unique = y_all(i_first_unique_sort);
% 
% arrow_dy = -y_last_unique + y_first_unique;
% arrow_du = -u_last_unique + u_first_unique;
% arraws = [u_last_unique',y_last_unique',arrow_du',arrow_dy'];

% 
% figure; hold on;
% plot(u_all,y_all,'-','color',[0.7,0.7,0.7])
% plot(0,0,'.r','markersize',6);
% plot(u_all,y_all,'b.','markersize',6)
% 
% for i=1:1:size(arraws,1)
%         headWidth = 5;
%         ah = annotation('arrow',...
%             'headStyle','cback1','HeadLength',5,'HeadWidth',headWidth,'linestyle','-');
%         set(ah,'parent',gca);
%         set(ah,'position',arraws(i,:));
% end
% ylim([-2.3,2.8])
% % for i = 1:1:c
% %     h1 = plot(mi(i,1),mi(i,2),'o','markerfacecolor',color(i,:));
% %     plot(x_sigma(:,i) + mi(i,1),y_sigma(:,i)+ mi(i,2), 'Color',color(i,:))
% %
% % end:-1:1) = p_all_unsorted(1:2);
% for i = 1:1:length(p_all_unsorted)
%     text(p_all_unsorted(i)+0.05,mi(i_closest(i),2),num2str(i),'color','r','FontSize',24)
% end
% %ylim([0,55])
% %set(h1, 'markerfacecolor', get(h1, 'color'));
% xlabel('Input $u$')
% ylabel('System output $y$')
% ax = gca;
% ax.Toolbar.Visible = 'off';
% set(ax,'fontname','Times', 'FontSize', 14);
% %h = title('Zaporedne točke vzbujanja');
% set(ax,'FontSize',14,'FontWeight','Normal','fontname','Times')
% h = legend('Consecutive excitation','Operating point','Measurement',...
%   'location','northwest');
% %h = legend('Excitation order','Sample','location','northwest');
% set(ax,'FontSize',14)
% name = ['FRLS_heuristic_zaporedne_tocke_povezane','.pdf'];
% exportgraphics(gcf, name,'ContentType','vector');

% figure;
% subplot(2,1,1); hold on;
% plot(phi_P_phi_posteriori(1:end-1),'b','linewidth',1.2)
% ylabel(' $\varphi P \varphi^T$')
% xlabel('Time step $k$')
% ax = gca;
% ax.Toolbar.Visible = 'off';
% set(ax,'fontname','Times', 'FontSize', 14);
% subplot(2,1,2);
% plot(alpha_f_all((1:end-1)),'b-.','linewidth',1.2)
% %title('Adaptation parameter')
% ylabel('Parameter $\alpha_f$')
% xlabel('Time step $k$')
% %xlim([1,length(phi_P_phi_posteriori)])
% ax = gca;
% ax.Toolbar.Visible = 'off';
% set(ax,'fontname','Times', 'FontSize', 14);
% name = ['FRLS_heuristika_phi_p_phi','.pdf'];
%exportgraphics(gcf, name,'ContentType','vector');

% figure;
% for i = 1:1:c
%     subplot(c,1,i)
%     plot(real(z_i_all(:,:,i)),imag(z_i_all(:,:,i)),'.','color',color(i,:))
%     xlim([-1,1])
% end
% 
% figure; hold on;
% xunit = 1 * cos(0:0.01:2*pi);
% yunit = 1 * sin(0:0.01:2*pi);
% plot(xunit,yunit,'k')
% for i = 1:1:c
%     plot(real(z_i_all(1,:,i)),imag(z_i_all(1,:,i)),'o','Linewidth',5,'MarkerSize',5,...
%     'MarkerEdgeColor',color(i,:),'MarkerFaceColor',color(i,:))
%     plot(real(z_i_all(end,:,i)),imag(z_i_all(end,:,i)),'s','Linewidth',5,'MarkerSize',5,...
%     'MarkerEdgeColor',color(i,:),'MarkerFaceColor',color(i,:))
%     plot(real(z_i_all(:,:,i)),imag(z_i_all(:,:,i)),'.','Linewidth',10,'MarkerSize',10,...
%     'MarkerEdgeColor',color(i,:),'MarkerFaceColor',color(i,:))
% end
% 
% 
% figure; hold on;
% plot(w); axis manual
% plot(y_m_w)
% plot(y_m_wo)
% plot(z(:,2))
% xlabel('Time step')
% ylabel('Output signal')
% legend('Reference','Model output with','Model output without','System output','Location','south')
% ax = gca;
% ax.Toolbar.Visible = 'off';
% set(ax,'fontname','Times', 'FontSize', 14);
% 
% 
% 
% figure
% subplot(2,1,1);hold on
% plot(w); axis manual
% plot(y_m_w,'k')
% plot(y_m_wo,'k--')
% plot(z(:,2))
% xlabel('Time step')
% ylabel('Output signal')
% legend('Reference','Model output with','Model output without','System output','Location','south')
% ax = gca;
% ax.Toolbar.Visible = 'off';
% set(ax,'fontname','Times', 'FontSize', 10);
% 
% subplot(2,1,2);hold on
% plot(u)
% %plot(e_p)
% xlabel('Time step')
% ylabel('System input')
% %legend('System input','Control error','Location','best')
% ax = gca;
% ax.Toolbar.Visible = 'off';
% set(ax,'fontname','Times', 'FontSize', 10);
% name = ['PWL_MPC_identification_','control','.pdf'];
% exportgraphics(gcf, name,'ContentType','vector');
% 
% figure;hold on;
% plot(a_r_all)
% plot(n_tau_all)
% plot(horiz_all)
% 
% legend('a_{r}','n_{tau1}','n_{tau2}','horiz','')
% 

% figure; hold on;
% plot(horiz_all)
% ylabel('Time horizon')
% xlabel('Time step')
% ax = gca;
% ax.Toolbar.Visible = 'off';
% set(ax,'fontname','Times', 'FontSize', 14);
% name = ['PWL_MPC_identification_','horizon','.pdf'];
% 
% exportgraphics(gcf, name,'ContentType','vector');
% %load PWL_FRLS_MPC_identification

figure; hold on;
plot(z(:,2),'b')
plot(y_m_wo,'r')

