nKMCSteps = 10000;
K = load('K');
[vN, IN, P, vR] = rmat2kmc(K)
ind = 1;
indLast = ind;
NTrans = zeros(4); % record transitions between states 1..4
tTot = 0;
for i=1:nKMCSteps
  indOld = ind;
  [ind,t] = kmc1(ind, vN, IN, P, K, vR);
  tTot = tTot + t;
  if ismember(ind,[1:4])
    NTrans(ind,indLast) = NTrans(indLast,ind) +1;
    indLast = ind;
  end
end
save NTrans NTrans
