function vAssignment = regroupK(K,nClust)
% groups states with rate matrix K into nClust clusters
n=size(K,1);
vAssignment=[1:n]';
nNow = n;
Kmin = min(K,K');
Ujo = eye(n);
while nNow>nClust
  vkmin = max(Kmin);
  [ujo,i1] = max(vkmin);
  vkmin(i1) = 0;
  [ujo,i2] = max(vkmin);
  Ujo += Kmin>=ujo;
  Kmin -= Ujo*ujo;
  vi1 = find(vAssignment==vAssignment(i1));
  vi2 = find(vAssignment==vAssignment(i2));
  numMin = size(vi1,1)+size(vi2,1);
  SubMat = Ujo([vi1;vi2],[vi1;vi2]);
  if (sum(sum(SubMat))==numMin^2)
    vAssignment(vi2) = vAssignment(i1);
  endif
  nNow = size(unique(vAssignment),1);
endwhile
% make group numbers between 1 and the number of groups
vOld = vAssignment;
vUniq = unique(vAssignment);
for iClust=1:nClust
  vAssignment(find(vOld==vUniq(iClust))) = iClust;
endfor
