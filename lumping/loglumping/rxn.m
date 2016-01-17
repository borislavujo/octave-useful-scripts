function [lkAB, lkBA] = rxn(Kl,vlp,va,howFine,lstartdt,thresh)
% calculates rate constant between 2 states from a 3x3 rate matrix K and log of equilibrium populations vlp
% grouped states are states 1 and 2, rate between (1+2) and (3) is calculated
n = size(Kl,1); vlp = vlp - logSumExp(vlp);
Kl(find(Kl==0)) = -1e9; % zero rates
Kl([1:n+1:n^2]) = -1e9; % change diagonal rates to 9
if (nargin<4) howFine = 2; endif % each exponential step is divided into 2^howFine linear propagation steps
if (nargin<5) lstartdt = -6; endif % the initial time step is exp(lstartdt) x (time of the fastest process in the network)
if (nargin<6) thresh = -5; endif % criterion for stopping the propagation: 1-ratioEnd relaxed for all groupings
nFine = 2^howFine;
lds = -max(max(Kl))+lstartdt-howFine*log(2); % max dt calculated from greatest rate
Ls = rm2tm(Kl,lds);
[Ds,Ns] = tm2dn(Ls,vlp);
L = Ls; D = Ds; N = Ns;
for iFine=1:howFine [L,D,N] = multLogMat(vlp,L,D,N); endfor
X = [-1e9,0]; ldt = lds + howFine*log(2);
while (max(max(D))>thresh)
  x = lpopul(L,D,N,vlp,va); 
  X = [X; ldt, x];
  Lp = L; Dp = D; Np = N; ltnow = ldt;
  for iFine=1:nFine
    [Lp,Dp,Np]=multLogMat(vlp,Lp,Dp,Np,Ls,Ds,Ns);
    x = lpopul(Lp,Dp,Np,vlp,va);
    ltnow = logSumExp([ltnow;lds]);
    X = [X; ltnow, x];
  endfor
  [L,D,N] = multLogMat(vlp,L,D,N); % transition matrix for dt -> 2 dt
  [Ls,Ds,Ns] = multLogMat(vlp,Ls,Ds,Ns); % transition matrix for dt / 2^howFine -> dt / 2^(howFine-1)
  ldt += log(2); % double dt
  lds += log(2);
endwhile
ltau = intLog(X);
vypoctov = size(X,1)
pl1 = logSumExp(vlp(va==1));
lkBA = pl1-ltau;
lkAB = subFrom1(pl1)-ltau;
%ltau = expTail(vT,[Kl(1,3),Kl(2,3),Kl(3,1),Kl(3,2)],vlp,ltau); % add the tail (since the propagation ends early; see ratioEnd)
%[lkTSTAB,lkTSTBA] = tstlog(Kl,vlp); % calculate TST rates
%lkAB = logSumExp(vlp(1:2))-ltau;
%tstrxn = exp(lkTSTAB-lkAB)
%if (lkAB>lkTSTAB)
%  lkAB = lkTSTAB;
%  tstwatsmaller = 1
%endif
%lkBA = lkAB+vlp(3)-logSumExp(vlp(1:2));
endfunction

function ltau = intLog(X)
% log version of trapz
n = size(X,1);
ltau = -1e9;
for i=1:n-1
  ldt = logDiffExp(X(i+1,1),X(i,1));
  pl = logSumExp(X(i:i+1,2))-log(2);
  ltau = logSumExp([ltau,ldt+pl]);
endfor
endfunction