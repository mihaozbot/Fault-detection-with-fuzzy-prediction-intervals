
%Input-output sample for clustering
z(k,:) = [u(k),y(k), d(k)];

%Filtering
d_f(k) = sum(A_f)*z(k,3) - A_f(2:(n+1))*d_f(k-1:-1:k-n);
u_f(k) = sum(A_f)*z(k,1) - A_f(2:(n+1))*u_f(k-1:-1:k-n);
y_f(k) = sum(A_f)*z(k,2) - A_f(2:(n+1))*y_f(k-1:-1:k-n);

%Regression vektor(RLS)
if enable_disturbance
    phi = [u_f((k - 1):-1:(k - m))',-y_f((k - 1):-1:(k - n))', d_f(k-1), 1]';
else
    phi = [u_f((k - 1):-1:(k - m))',-y_f((k - 1):-1:(k - n))', 1]';
end