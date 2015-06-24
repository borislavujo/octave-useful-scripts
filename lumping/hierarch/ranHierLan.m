function [K, vPop] = ranHierLan(nLev, beta, scaleTS, hierParam)
% creates a random hierarchical landscape and stores it in terms of a kinetic network
nMins = 2^nLev;
vPop = rand(nMins,1); vPop = vPop/sum(vPop); % populations - uniform random numbers
vMinF = -log(vPop)/beta; % free energies of minima
% probability distribution of relaxation times is exp(scaleTS*x), where x=rand(1)*deltaLev (deltaLev is the number of dividing layers)
Trxn = zeros(nMins);
K = zeros(nMins);
vdi = 2.^[0:nLev-1]';
for i1=1:nMins-1
  for i2=i1+1:nMins
    logtau = scaleTS*(rand(1)^hierParam)*sum(ceil(i1./vdi)~=ceil(i2./vdi));
    K(i2,i1) = (vPop(i2)/(vPop(i1)+vPop(i2)))/exp(logtau);
    K(i1,i2) = (vPop(i1)/(vPop(i1)+vPop(i2)))/exp(logtau);
  endfor
endfor
K = K - diag(sum(K,1));
