% Assuming 'c' is the total number of clusters
for jo = 1:c

    % Recursively identify the parameters of the consequence
    gamma = (P(:,:,jo)*phi)/(phi'*P(:,:,jo)*phi + 1/NGamma(jo));
    P(:,:,jo) = (eye(n_phi) - gamma*phi')*P(:,:,jo);
    e_RLS(k) = y_f(k) - phi'*theta(:,jo);
    theta(:,jo) = theta(:,jo) + gamma*e_RLS(k);

    e_phi = [phi' - phi_c(jo,:)]'; % Distance of new data from center
    phi_c(jo,:) = phi_c(jo,:) + NGamma(jo)/(1 + n_c(jo))*e_phi'; % Center update

    %     % Compute projection
    %     Sigma(:,:,jo) = S_c(:,:,jo)/n_c(jo);
    %     [~,S,V] = svd(Sigma(:,:,jo));
    %     delta(:,:,ji) = diag([max(abs(V(1,:)'.*diag(sqrt(S)))),...
    %         max(abs(V(2,:)'.*diag(sqrt(S))))]);

    % % Incremental clustering of the antecedent clusters
    e_c = [z(k,1:n_z) - mi(jo,:)]'; % Distance of new data from center
    mi(jo,:) = mi(jo,:) + NGamma(jo)/(1 + n_c(jo))*e_c'; % Center update
    S_c(:,:,jo) = S_c(:,:,jo) + NGamma(jo)*e_c*(z(k,1:n_z) - mi(jo,:)); % Un-normalized covariance matrix
    n_c(jo) = n_c(jo) + NGamma(jo); % Increase number of samples in cluster

end

% % % Incremental clustering of the antecedent clusters
% e_c = [z(k,1:n_z) - mi(jg,:)]'; % Distance of new data from center
% mi(jg,:) = mi(jg,:) + 1/(1 + n_c(jg))*e_c'; % Center update
% S_c(:,:,jg) = S_c(:,:,jg) + e_c*(z(k,1:n_z) - mi(jg,:)); % Un-normalized covariance matrix
% n_c(jg) = n_c(jg) + 1; % Increase number of samples in cluster
