function [Nxy,Pxy,vpxx,vnpxx,vtaux,vn] = rm2ngt(Kl)
% transforms a log rate matrix into structures required by a 
% Nxy - indices of neighbours
% Pxy - log probabilities to neighbours
% vpxx - log probabilities to self
% vtaux - log waiting times
% vn - numbers of neighbours
malo = -1e9; nula = -99e9;
n = size(Kl,1)
vpxx = nula*ones(n,1);
vnpxx = zeros(n,1);
vn = zeros(n,1);
for i=1:n
  vn(i) = size(find(Kl(i,:)>malo),2);
endfor
mn = max(vn);
vtaux = -nula*ones(n,1);
Nxy = zeros(mn,n);
Pxy = nula*ones(mn,n);
for i=1:n
  if (vn(i)>0)
    vnnow = find(Kl(:,i)>malo)';
    Nxy(1:vn(i),i) = vnnow;
    vknow = Kl(vnnow,i)
    tauxnow = -logSumExp(vknow);
    vtaux(i) = tauxnow;
    vpnow = vknow + vtaux(i);
    Pxy(1:vn(i),i) = vpnow;
  else
    vtaux(i) = -nula; % -nula = infty
  endif
endfor
