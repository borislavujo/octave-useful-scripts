function [vNeighbours, IndNeighbours, Proportions, vRates] = rmat2kmc(K)
% makes matrices useful for KMC from a rate matrix
n = size(K,1);
vNeighbours = zeros(n,1);
IndNeighbours = sparse(n-1,n);
Proportions = sparse(n-1,n);
K = K - diag(diag(K));
for i=1:n
  vNi = find(K(:,i));
  nN = size(vNi,1);
  vNeighbours(i) = nN;
  IndNeighbours(1:nN,i) = vNi;
  scale = sum(K(:,i));
  s = 0;
  for j=1:nN
    s = s + K(IndNeighbours(j,i),i);
    Proportions(j,i) = s/scale;
  end
end
vRates = sum(K);