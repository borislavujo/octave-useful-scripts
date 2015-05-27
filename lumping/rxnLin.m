function tau = rxnLin(K,vpe,vasi)
% calculates relaxation time
% K .. rate matrix
% vp .. vector of equilibrium populations
% vasi .. assignment into states
ratioEnd = 1e-3;
stepSize = 0.1;

% checks
vujo = unique(vasi);
if (size(vujo,1)~=2) error('this calculates relaxation between 2 states only'); endif
n = size(K,1);
p1 = sum(vpe(find(vasi==vujo(1))))
p2 = sum(vpe(find(vasi==vujo(2))))
if (abs(p1+p2-1)>1e-3) error('sum(pop) must be 1'); endif
p = 1; vp = vpe.*(vasi==vujo(1)); vp = vp / p;

% key variables for propagation
dt = stepSize/max(abs(diag(K))) % quite arbitrary compromise between efficiency and accuracy
F = (eye(n)+dt*K);
t = 0;
iStep = 0;
nCheck = 1;
tau = 0;
while (p>ratioEnd)
  vp = F*vp;
  vp = vp/sum(vp);
  pold = p;
  p = (sum(vp(find(vasi==vujo(1))))-p1)/(1-p1);
  iStep++;
  t+=dt;
  tau += dt * (pold+p)/2;
endwhile

tauNow = tau
dp = sum((K*vp)(find(vasi==vujo(1))))/(1-p1)
k = -dp/p
tau = tau + p/k
