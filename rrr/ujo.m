% test script for the RRR
[Kl, vpl] = rlm(15,2,2); % generate a random (full) log rate matrix
[vi,FEl,SSl,TSTl,NFl] = rrr(Kl,vpl,0,[ones(7,1);2*ones(8,1)]) % group into a 2-state model with arbitrary observables
[vi,FEl,SSl,TSTl,NFl] = rrr(Kl,vpl,2) % group into a 2-state model with observables defined by the grouping
ratioMaxMin = FEl(1,2)/SSl(1,2)
minTSTError = TSTl(1,2)/FEl(1,2)
ratioNF2TST = NFl(1,2)/TSTl(1,2)

