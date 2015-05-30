tic
ratioEnd = 1e-5;
startdt = 1e-1;
WK = load('K','-ascii'); % K is a sparse rate matrix
n = max(max(WK(:,1:2)))
K = zeros(n) ; % full rate matrix
for i=1:size(WK,1)
    i1 = WK(i,1);
    i2 = WK(i,2);
    k = WK(i,3);
    K(i1,i2) = k;
end
vpe = load('vp'); % equilibrium population
vasi = load('vas'); % assignment into states

p1 = sum(vpe(find(vasi==1))) % equilibrium population of state 1
p2 = sum(vpe(find(vasi==2)))
p = 1; vp0 = vpe.*(vasi==1); vp0 = vp0 / p;

dt0 = startdt/max(abs(diag(K))) % quite arbitrary compromise between efficiency and accuracy
F = (eye(n)+dt0*K); % transition matrix
t = 0; tau = 0; progress = 1-p
beforePropagation = toc % time this
while (1-progress>ratioEnd)
  vp = F*vp0; % calculate population at time t
  sparseness = nnz(F)/(n^2) % display transition matrix 'sparseness'
  F = F*F; % t -> 2*t
  vp = vp/sum(vp); % correct for numerical errors
% fit to exponential
  pold = p;
  p = (sum(vp(find(vasi==1)))-p1)/(1-p1)
  a0 = pold^2/p
  k = log(pold/p)/t
  tau = tau + (pold-p)/(k*a0) % integrate the segment
  tauEstd = tau + p/k % tau estimated from current tau and k
  progress = 1-p
  if (t==0) t=dt0; else t = t*2; end
  computerTime = toc
end

tauNow = tau
vdp = K*vp;
dp = sum(vdp(find(vasi==1)))/(1-p1)
k = -dp/p
tau = tau + p/k % add the tail
save -ascii tauFinal tau
