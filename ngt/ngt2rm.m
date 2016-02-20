function Kl = ngt2rm(Nxy,Pxy,vpxx,vnpxx,vtaux,vn)
% transforms a log rate matrix into structures required by a 
% Nxy - indices of neighbours
% Pxy - log probabilities to neighbours
% vpxx - log probabilities to self - not used
% vtaux - log waiting times
% vn - numbers of neighbours
n = size(vn,1)
Kl = -99e9*ones(n);
for ib = 1:n
  if (vn(ib)>0)
    for jxb = 1:vn(ib)
      ig = Nxy(jxb,ib);
      Kl(ig,ib) = Pxy(jxb,ib)-vtaux(ib);
    endfor
  endif
endfor
