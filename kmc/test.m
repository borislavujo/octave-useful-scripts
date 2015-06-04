%WK = load('K','-ascii'); % K is a sparse rate matrix
%n = max(max(WK(:,1:2)))
%K = zeros(n) ; % full rate matrix
%for i=1:size(WK,1)
%    i1 = WK(i,1);
%    i2 = WK(i,2);
%    k = WK(i,3);
%    K(i1,i2) = k;
%end
K = load('K');

[vN, IN, P] = rmat2kmc(K)
ind = 1;
indLast = ind;
NTrans = zeros(4);
tTot = 0;
for i=1:100000
  indOld = ind;
  [ind,t] = kmc1(ind, vN, IN, P, K)
  tTot = tTot + t
  if ismember(ind,[1:4])
    NTrans(ind,indLast) = NTrans(indLast,ind) +1
    indLast = ind;
  endif
end
save NTrans NTrans