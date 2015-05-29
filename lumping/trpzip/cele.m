Mins = load('min.data');
TS = load('ts.data');
nm = size(Mins,1)
vp = makePop(Mins);
K = makeK(Mins,TS);

vas1 = load('vas1')
vasi1 = 2*ones(nm,1);
vasi1(vas1) = 1;
tau1 = rxnLin(K,vp,vasi1)
save tau1 tau1

vas4 = load('vas4')
vasi4 = 2*ones(nm,1);
vasi4(vas4) = 1;
tau4 = rxnLin(K,vp,vasi4)
save tau4 tau4

vas6 = load('vas6')
vasi6 = 2*ones(nm,1);
vasi6(vas6) = 1;
tau6 = rxnLin(K,vp,vasi6)
save tau6 tau6

vas10 = load('vas10')
vasi10 = 2*ones(nm,1);
vasi10(vas10) = 1;
tau10 = rxnLin(K,vp,vasi10)
save tau10 tau10

