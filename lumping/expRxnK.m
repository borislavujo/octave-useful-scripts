function [kAB, kBA] = expRxnK(K,va,vpe, startdt, ratioEnd)
if (nargin<4) startdt = 1.5e-2; endif% the initial time step is startdt x (time of the fastest process in the network)
if (nargin<5) ratioEnd = 1e-3; endif % criterion for stopping the propagation: 1-ratioEnd relaxed for all groupings
n = size(K,1);
vpe = vpe/sum(vpe); P = diag(vpe); % normalise equilibrium populations
vpei = ones(size(vpe))./vpe; Pi = diag(vpei); % inverse equilibrium population
G = Pi*K; K = P*(G+G')/2; % symmetrise the rate matrix
K = K - diag(sum(K)); % just to ensure the rate matrix is correct
p1 = sum(vpe.*va); p2 = 1-p1; % equilibrium population of states 1 and 2
vp0 = vpe.*va/p1; p = (sum(vp0)-p1)./p2; % p is 1 at the beginning
dt0 = startdt/max(abs(diag(K))); % quite arbitrary compromise between efficiency and accuracy
G = (eye(n)+dt0*K); % transition matrix
t = dt0/2; tau = 0; 
% relaxation
while (p>ratioEnd)
  K = Pi*G; G = P*(K+K')/2; % 
  F = G + diag(ones(1,n)-sum(G,1)); % sum of elements in each column must be exactly 1
  F(abs(F) < sqrt(realmin)) = 0; % avoid complicated algebraic manipulations with numbers of order 1e-154 and less
  vp = F*vp0; % calculate population at time t
  pold = p; p = (sum(va.*vp)-p1)./p2; 
  tau = tau + t*(p+pold)/2;
  t = t*2; % t -> 2*t
  G = F*F; % t -> 2*t
end

vdp = K*vp; dp = sum(vdp.*va)./p2; k = -dp./p; 
tau = (tau + p./k); kAB = p2/tau; kBA = p1/tau; % add the tail
