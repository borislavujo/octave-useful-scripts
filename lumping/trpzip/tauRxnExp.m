% difference from the others (g1,g4,g6)
% sparse matrix is used instead of a full matrix

% calculates relaxation time
% K .. rate matrix
% vp .. vector of equilibrium populations
% vasi .. assignment into 2 states
% vConn .. indices of states connected to the native state
tic
ratioEnd = 1e-5;
startdt = 1e-1;
%F = load('F');
K = load('K').K;
%K = K/max(max(K));
%save K K
maxTerm = max(max(K))
minTerm = min(min(K))
%K = full(K);
n = size(K,1)
% remove disconnected states
vConn = load('vConn');
K = K(vConn,:);
K = K(:,vConn);
[n,m] = size(K)
vpe = load('vp');
vpe = vpe(vConn);
vpe = vpe/sum(vpe);
% assignment to groups
vas = load('vas');
vasi = 2*ones(n,1);
for i=1:n
  if (ismember(vConn(i),vas))
    vasi(i) = 1;
  endif
endfor

% checks
p1 = sum(vpe(find(vasi==1)))
p2 = sum(vpe(find(vasi==2)))
p = 1; vp0 = vpe.*(vasi==1); vp0 = vp0 / p;

% key variables for propagation
dt0 = startdt/max(abs(diag(K))) % quite arbitrary compromise between efficiency and accuracy
F = (eye(n)+dt0*K);
t = 0;
tau = 0;
progress = 1-p
beforePropagation = toc
while (1-progress>ratioEnd)
  vp = F*vp0;
  sparseness = nnz(F)/(n^2)
  F = F*F;
  vp = vp/sum(vp);
  pold = p;
  p = (sum(vp(find(vasi==1)))-p1)/(1-p1)
  a0 = pold^2/p
  k = log(pold/p)/t
  tau += (pold-p)/(k*a0)
  tauEstd = tau + p/k
  progress = 1-p
  if (t==0) t=dt0 else t *= 2 endif
  computerTime = toc
endwhile

tauNow = tau
dp = sum((K*vp)(find(vasi==1)))/(1-p1)
k = -dp/p
tau = tau + p/k
