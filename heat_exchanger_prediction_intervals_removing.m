
%Remove check, so that multiple clusters are removed at the end
remove_check = zeros(1,c) > 0;

% %Remove unstable clusters
% if any(abs(z_i)>(1-1e-4))
%      remove_check(c) = 1;
%      disp('Remove unstable cluster')
% end

%Removing based on overlapping, not applicable here


%Removal of rules based on their MSE
remove_check((MSE > kappa_remove)) = 1;

%Inforamtion on the removal of rules
if sum(remove_check) > 1
    disp('Multiple removal')
elseif any(remove_check)
    disp(['Remove ',num2str(find(remove_check))])
end

%Removal mechanism
if any(remove_check)

    n_c(remove_check) = [];
    mi(remove_check,:) = [];
    S_c(:,:,remove_check) = [];
    theta(:,remove_check) = [];
    b_fuzzy(remove_check) = [];
    c = c-sum(remove_check);
    MSE(remove_check) = [];
    P(:,:,remove_check) = [];
    %Sigma(:,:,remove_check) = [];
    %delta(:,:,remove_check) = [];
    Gamma(remove_check) = [];
    phi_c(c,:) = [];
end
