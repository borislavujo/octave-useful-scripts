function [Kl, ltau] = allngt(Kl,vi)
if (size(Kl,1)>2)
  [Nxy,Pxy,vpxx,vnpxx,vtaux,vn] = rm2ngt(Kl);
  n = size(Kl,1);
  vred = setdiff([1:n]',vi);
  nred = size(vred,1)
  while (nred>0)
% start from the ones with minimum number of contacts
    [kelo,kery] = min(vn(vred))
    ired = vred(kery)
    [Nxy,Pxy,vpxx,vnpxx,vtaux,vn] = ngt1(Nxy,Pxy,vpxx,vnpxx,vtaux,vn,ired);
    vred(kery) = [];
% reindex vred
    if and(size(vred,1)>0,size(vred,2)>0)
      vysoke = find(vred>kery);
      if and(size(vysoke,1)>0,size(vysoke,2)>0)
	for i=vysoke'
	  vred(i) = vred(i)-1;
	endfor
      endif
    endif
    nred--;
  endwhile
  Kl = ngt2rm(Nxy,Pxy,vpxx,vnpxx,vtaux,vn);
endif
ltau = -logSumExp([Kl(1,2);Kl(2,1)]);
endfunction