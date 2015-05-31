tic
ratioEnd = 1e-5;
startdt = 1e-1;
% load the rate matrix
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
Lump = load('Lump','-ascii'); % assignment into states
nLumps = size(Lump,2);
% set up the starting population
Peq = repmat(vpe,1,nLumps);
vp1 = sum(Peq.*Lump) % equilibrium population of state 1
Pop0 = Peq.*Lump./repmat(vp1,n,1);
vp = sum(Pop0);
% the first transition matrix
dt0 = startdt/max(abs(diag(K))) % quite arbitrary compromise between efficiency and accuracy
F = (eye(n)+dt0*K); % transition matrix
t = 0; vTau = 0; beforePropagation = toc % time this
% relaxation
while (max(abs(vp))>ratioEnd)
  Pop = F*Pop0; % calculate population at time t
  Pop = Pop ./ repmat(sum(Pop),n,1);
  F = F*F; % t -> 2*t
  F(abs(F) < sqrt(realmin)) = 0;
  sparseness = nnz(F)/(n^2) % display transition matrix 'sparseness'
% fit to exponential
  vpold = vp;
  p = sum(Lump.*Pop)./(ones(1,nLump)-vp1)
  vTau = vTau + t*(vp*vpold)/2
% estimate final tau
  DP = K*Pop;
  vdp = sum(DP.*Lump)./(ones(1,nLump)-vp1)
  vk = -vdp./vp
  vTauEstd = vTau + vp./vk % tau estimated from current tau and k
  if (t==0) t=dt0; else t = t*2; end
  computerTime = toc
end

DP = K*Pop;
vdp = sum(DP.*Lump)./(ones(1,nLump)-vp1)
vk = -vdp./vp
vTau = vTau + vp./vk % add the tail
save -ascii vTau vTau
