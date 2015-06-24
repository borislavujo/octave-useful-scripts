tic
time0 = time;
ratioEnd = 1e-3; % criterion for stopping the propagation: 1-ratioEnd relaxed for all groupings
startdt = 1e-2; % the initial time step is startdt x (time of the fastest process in the network)
scalek = 1; % to mitigate numerical errors
odkal = 1; % select which groupings you are interested in
pokal = 20; % select which groupings you are interested in
nHow = 4;
% load the rate matrix .. compatible with matlab
WK = load('K','-ascii'); % K is a sparse rate matrix
n = max(max(WK(:,1:2)))
K = zeros(n) ; % full rate matrix
for i=1:size(WK,1)
    i1 = WK(i,1); i2 = WK(i,2); k = WK(i,3);
    K(i1,i2) = k/scalek;
end
vpe = load('vp'); % equilibrium population
Lump = load('Lump','-ascii')(:,odkal:pokal); % assignment into states
nLumps = size(Lump,2);
vkba = zeros(1,nLumps);
for iLump=1:nLumps
  [kAB, kBA] = rxnK(K,Lump(:,iLump),vpe,nHow,startdt,ratioEnd)
  vkba(iLump) = kBA;
endfor

Kngt = load('ktst');
hold on; 
title ("effect of regrouping on folding rate"); 
xlabel("regrouping threshold / [kcal/mol]"); 
ylabel("log10(folding rate / [1/s])"); 
plot(Kngt(odkal:pokal,1),log(vkba)/log(10),"-;relaxation-based lumping;")
plot(Kngt(odkal:pokal,1),log(Kngt(odkal:pokal,2))/log(10),"-.;TST-based lumping;"); 
hold off;
print -depsc regroup.eps

