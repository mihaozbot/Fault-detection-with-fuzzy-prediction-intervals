%The parameter vector With and WithOut the last added cluster
%Theta_w = theta(:, c);
Theta_wo = theta(:, jg);
icg = setdiff(1:c, jg);% Create an index vector that excludes 'jg'
if c > 1
    % Use these indices to compute the weighted sum
    Theta_wo = sum((repmat(Gamma(icg)', n_phi, 1) .* theta(:, icg)) / sum(Gamma(icg)), 2);
end

%Construct the regression vector at time step k
if enable_disturbance
    %Phi_m_w = [z(k-1:-1:k-m,1)',-y_m_w(k-1:-1:k-n)',z(k-1,3),1]';
    Phi_m_wo = [z(k-1:-1:k-m,1)',-y_m_wo(k-1:-1:k-n)',z(k-1,3),1]';
else
    %Phi_m_w = [z(k-1:-1:k-m,1)',-y_m_w(k-1:-1:k-n)',1]';
    Phi_m_wo = [z(k-1:-1:k-m,1)',-y_m_wo(k-1:-1:k-n)',1]';
end
%y_m_w(k) = Phi_m_w'*Theta_w;
y_m_wo(k) = Phi_m_wo'*Theta_wo;
if isnan(y_m_wo(k))
    y_m_wo(k) = y_m_wo(k-1);
end

%Compute the predicted model output at time step k
Theta = Theta_wo;
y_m_p = y_m_wo(k)';

%Compute the model error at time step k
%e_est_w(k) = z(k,2)-y_m_w(k);
e_est_wo(k) = z(k,2)-y_m_wo(k);

%Recursivelly compute the mean square error of the error
%Compute the fuzzy model error
d_MSE_wo = (e_est_wo(k)^2)*(Gamma(1:c)./sum(Gamma(icg))).^2;
for g = 1:c
    if g==jg
        continue
    end

    %Update the fuzzy mean squared error (fMSE)
    if b_fuzzy(g,1)>(n_phi+1)
        MSE(g,1) = (1./(b_fuzzy(g,1) + (...
            Gamma(g,1)/sum(Gamma(icg,1)))^2 - (n_phi+1)))*...
            ((b_fuzzy(g,1) - (n_phi+1))*MSE(g,1) + d_MSE_wo(g,1));

    end

    %Update the fuzzy acumulated membership
    b_fuzzy(g,1) = b_fuzzy(g,1) + (Gamma(g,1)/sum(Gamma(icg,1))).^2;

end

%Debugging checks ...
if sum(b_fuzzy) > 2*k
    disp('ERROR! b_fuzzy can not be this large!')
    disp(b_fuzzy)
end
if any(MSE<0)
    disp('ERROR! MSE is negative')
    disp(MSE)
end
if any(isnan(MSE))
    disp('ERROR! MSE is NaN')
    disp(MSE)
end


