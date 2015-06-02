tic
time0 = time;
ratioEnd = 1e-3; % criterion for stopping the propagation: 1-ratioEnd relaxed for all groupings
startdt = 1e-2; % the initial time step is startdt x (time of the fastest process in the network)
scalek = 1000; % to mitigate numerical errors
odkal = 1; % select which groupings you are interested in
pokal = 20; % select which groupings you are interested in
% load the rate matrix .. compatible with matlab
WK = load('K','-ascii'); % K is a sparse rate matrix
n = max(max(WK(:,1:2)))
K = zeros(n) ; % full rate matrix
for i=1:size(WK,1)
    i1 = WK(i,1); i2 = WK(i,2); k = WK(i,3);
    K(i1,i2) = k/scalek;
end
vpe = load('vp'); % equilibrium population
vpei = ones(size(vpe))./vpe; % inverse equilibrium population
Lump = load('Lump','-ascii')(:,odkal:pokal); % assignment into states
nLumps = size(Lump,2);
% symmetrisation of the matrix
K = K - diag(sum(K)); % just to ensure the rate matrix is correct, so all rates in a column sum to 0
G = diag(vpei)*K; % G is just a temporary storage; a symmetric matrix
F = min(abs(G),abs(G'))/2; % F is used as a temporary array here; minimum instead of (G+G')/2 for better numerical stability
% this symmetrisation is a bit suspicious
K = diag(vpe)*F;
K = K - diag(sum(K)); % just to ensure the rate matrix is correct
% set up the starting population
Peq = repmat(vpe,1,nLumps);
vp1 = sum(Peq.*Lump) % equilibrium population of state 1
vp2 = ones(1,nLumps)-vp1;
Pop0 = Peq.*Lump./repmat(vp1,n,1);
vp = (sum(Pop0)-vp1)./vp2
% the first transition matrix
dt0 = startdt/max(abs(diag(K))) % quite arbitrary compromise between efficiency and accuracy
G = (eye(n)+dt0*K); % transition matrix
testMustBe0 = min(min(G))
testMustBeLT1 = max(max(G))
t = 0; vTau = 0; beforePropagation = toc % time this
% relaxation
while (max(abs(vp))>ratioEnd)
% numerical corrections
  F = diag(vpei)*G; 
  G = diag(vpe)*min(abs(F),abs(F')); % inv(P)*T = T'*inv(P) - instead of F=(F+F')/2, which is numerically instable
  testMustBe0 = min(min(G))
  testMustBeLT1 = max(max(G))
  [kolko,ind]=find(G>1)
  [ujo1,ujo2]=ind2sub([n,n],ind)
  vCorrect = ones(1,n)-sum(G,1);
  vpici = sort(vCorrect);
  vpici = [vpici(1:10);vpici(n-9:n)]
  totalCorrection = sum(abs(vCorrect))
  F = G + diag(vCorrect); % sum of elements in each column must be exactly 1
  F(abs(F) < sqrt(realmin)) = 0; % avoid complicated algebraic manipulations with numbers of order 1e-154 and less
  sparseness = nnz(F)/(n^2) % display transition matrix 'sparseness'
  Pop = F*Pop0; % calculate population at time t
  tt0 = time-time0
  G = F*F; % t -> 2*t
  timeOfMultiplic = time - time0 - tt0 % computer time
% fit to exponential
  vpold = vp;
  vp = (sum(Lump.*Pop)-vp1)./vp2
  vTau = vTau + t*(vp+vpold)/2
% estimate final tau
  DP = K*Pop;
  vdp = sum(DP.*Lump)./vp2
  vk = -vdp./vp
  vTauEstd = vTau + vp./vk % tau estimated from current tau and k
  if (t==0) t=dt0; else t = t*2 end
  computerTime = toc
end

DP = K*Pop;
vdp = sum(DP.*Lump)./vp2
vk = -vdp./vp
vTau = (vTau + vp./vk)/scalek % add the tail
save -ascii vTau vTau

vkba = vp1./vTau % folding rate
Kngt = load('ktst');
hold on; 
title ("effect of regrouping on folding rate"); 
xlabel("regrouping threshold / [kcal/mol]"); 
ylabel("log10(folding rate / [1/s])"); 
plot(Kngt(odkal:pokal,1),log(vkba)/log(10),"-;relaxation-based lumping;")
plot(Kngt(odkal:pokal,1),log(Kngt(odkal:pokal,2))/log(10),"-.;TST-based lumping;"); 
hold off;
print -depsc regroup.eps

