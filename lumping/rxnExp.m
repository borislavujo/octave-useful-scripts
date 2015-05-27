function tau = rxnExp(K,vpe,vasi)
% calculates relaxation time
% K .. rate matrix
% vp .. vector of equilibrium populations
% vasi .. assignment into 2 states
ratioEnd = 1e-5;
startdt = 1e-3;

% checks
vujo = unique(vasi);
if (size(vujo,1)~=2) error('this calculates relaxation between 2 states only'); endif
n = size(K,1);
p1 = sum(vpe(find(vasi==vujo(1))))
p2 = sum(vpe(find(vasi==vujo(2))))
if (abs(p1+p2-1)>1e-3) error('sum(pop) must be 1'); endif
p = 1; vp0 = vpe.*(vasi==vujo(1)); vp0 = vp0 / p;

% key variables for propagation
dt0 = startdt/max(abs(diag(K))) % quite arbitrary compromise between efficiency and accuracy
F = (eye(n)+dt0*K);
t = 0;
tau = 0;
progress = 1-p;
while (1-progress>ratioEnd)
  vp = F*vp0;
  F = F*F;
  vp = vp/sum(vp);
  pold = p;
  p = (sum(vp(find(vasi==vujo(1))))-p1)/(1-p1);
  a0 = pold^2/p;
  k = log(pold/p)/t;
  tau += (pold-p)/(k*a0);
  progress = 1-p;
  if (t==0) t=dt0; else t *= 2; endif
endwhile

tauNow = tau
dp = sum((K*vp)(find(vasi==vujo(1))))/(1-p1)
k = -dp/p
tau = tau + p/k
