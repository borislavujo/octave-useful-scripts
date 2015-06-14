nLev = 4;
n = 2^nLev; 
beta = 1; % temperature
scaleTS = 4;
hierParam = 0.9; % low value means that landscape is very hierarchical
nl = 2;
nTests = 1000;
sTests = 0;

for iTest=1:nTests
[K, vPop] = ranHierLan(nLev, beta, scaleTS, hierParam);
K1 = K(1:n/2,1:n/2);
K2 = K(n/2+1:n,n/2+1:n);
v1 = vPop(1:n/2);
v2 = vPop(n/2+1:n);
va1 = regroupK(K1,nl);
va2 = regroupK(K2,nl);
va2 = va2 + nl;
va = [va1; va2];
% calculate contact matrix
CM = zeros(2*nl);
for i1=1:2*nl
  for i2=1:2*nl
    if (nnz(K(find(va==i1),find(va==i2)))>0)
      CM(i1,i2) = 1;
    endif
  endfor
endfor

nHow = 3;
% calculate rate constants for each pair of states
vpsmall = zeros(2*nl,1);
for i=1:2*nl vpsmall(i) = sum(vPop(va==i)); endfor
vpsmall = vpsmall / sum(vpsmall);
Ksmall = zeros(2*nl);
for i1=1:2*nl-1
  for i2=i1+1:2*nl
% rate matrix is a block matrix
    vi1 = find(va==i1);
    vi2 = find(va==i2);
    nn1 = size(vi1,1);
    nn2 = size(vi2,1);
    vasi = [ones(nn1,1); zeros(nn2,1)];
    vpnow = [vPop(vi1); vPop(vi2)];
    Know = [K(vi1,vi1),K(vi1,vi2);K(vi2,vi1),K(vi2,vi2)];
    Know = Know - diag(sum(Know));
    [kab, kba] = rxnK(Know,vasi,vpnow,nHow);
    Ksmall(i1,i2) = kba;
    Ksmall(i2,i1) = kab;
  endfor
endfor
Ksmall = Ksmall - diag(sum(Ksmall));

% hierarchical
[kabsmall, kbasmall] = rxnK(Ksmall,[ones(nl,1);zeros(nl,1)],vpsmall,nHow)
% at once
[kabbig, kbabig] = rxnK(K, [ones(n/2,1);zeros(n/2,1)],vPop,nHow)
testMustBe1 = (kabsmall/kbasmall)/(kabbig/kbabig);
ratioFriction = kabbig/kabsmall
sTests += ratioFriction;
averRatio = sTests/iTest
endfor