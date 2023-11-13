function olap = bhattacharyya(S_c, n_c, mi, c, ic)

n_z = size(mi,2);

Sigma = reshape(cell2mat(...
    arrayfun(@(i) S_c(:,:,i)./n_c(i), 1:c,'UniformOutput', false)),...
    [n_z,n_z,c]);

det_inv_cov_i = arrayfun(@(i) det(Sigma(:,:,i)), 1:c);

olap = Inf(1,c);
for i = ic:1:c

    cov = (Sigma(:,:,ic)+Sigma(:,:,i))/2;
    olap(i) = (1/8)*(mi(ic,:) - mi(i,:))*inv(cov)*...
        (mi(ic,:) - mi(i,:))'...
        + (1/2)*log(det(cov)/(sqrt(det_inv_cov_i(ic)*det_inv_cov_i(i))));

end

olap(ic) = Inf;
olap(olap<0) = Inf;

if any(olap<0)
    disp('Overlap error olap < 0')
end
end

