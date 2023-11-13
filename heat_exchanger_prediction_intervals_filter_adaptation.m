%Compute the confidence interval
phi_P_phi(jg) = phi'*P(:,:,jg)*phi;

%Process gain and roots
A_p = [1,theta(m+1:m+n, jg)']';
B_P = [theta(1:m, jg)']';
K(k) = sum(B_P)/sum(A_p);
z_i = roots(A_p'); %Compute roots, should be improved tbh.

%Check if model is stable and has low confidence interval
if (max(abs(z_i))<1) && (phi_P_phi(jg)<kappa_c)
    %A_f = A_p; %Adapt the filter
end



