function [Nxy,Pxy,vpxx,vnpxx,vtaux,vn] = ngt1(Nxy,Pxy,vpxx,vnpxx,vtaux,vn,ix)
% function [Nxyn,Pxyn,vpxxn,vtauxn,vnn] = ngt1(Nxy,Pxy,vpxx,vtaux,vn,ix)
% reduces a network by state ired
% according to Wales 2009, JCP 130, 204111
% Nxy - indices of neighbours
% Pxy - log probabilities to neighbours
% vpxx - log probabilities to self
% vnpxx - log (1 - probabilities to self)
% vtaux - log waiting times
% vn - numbers of neighbours
%Nxyn = Nxy; Pxyn = Pxy; vpxxn = vpxx; vtauxn = vtaux; vnn = vn;
vnow = Nxy(1:vn(ix),ix);
nnow = size(vnow,1);
% 2 types of indices - real (stars w i) and in p (start w j)
% for each elimination, there are 3 states - x, b, g
% Nxy(ix,jxb) = ib etc.
if (nnow>0)
  for jbx=1:nnow
    ib = vnow(jbx)
    jxb = find(Nxy(:,ib)==ix);
    if (vpxx(ix)<-0.5)
      vtaux(ib) = logSumExp([vtaux(ib);Pxy(jxb,ib)+vtaux(ix)-logDiffExp(0,vpxx(ix))]); % (11) in Wales2009
    else
      vtaux(ib) = logSumExp([vtaux(ib);Pxy(jxb,ib)+vtaux(ix)-vnpxx(ix)]); % (11) in Wales2009
    endif
    for jgx=1:nnow
      ig = vnow(jgx)
      if (jbx==jgx)
	vpxx(ib) = logSumExp([vpxx(ib);Pxy(jbx,ix)+Pxy(jxb,ib)-logDiffExp(0,vpxx(ix))]); % (8) in Wales2009, for b=g
	if (vpxx(ix)<-0.5)
	  vnpxx(ib) = subFrom1(vpxx(ib));
	else
	  vnpxx(ib) = logDiffExp(vnpxx(ib),Pxy(jbx,ix)+Pxy(jxb,ib)-vnpxx(ix)); % (8) in Wales2009, for b=g
	  vpxx(ib) = subFrom1(vnpxx(ib));
	endif
      elseif ismember(ib,Nxy(1:vn(ig),ig))
	jxg = find(Nxy(:,ig)==ix);
	jbg = find(Nxy(:,ig)==ib);
	jgb = find(Nxy(:,ib)==ig);
	if (vpxx(ix)<-0.5)
	  Pxy(jgb,ib) = logSumExp([Pxy(jgb,ib);Pxy(jgx,ix)+Pxy(jxb,ib)-logDiffExp(0,vpxx(ix))]); % (8) in Wales2009, for b~=g
	else
	  Pxy(jgb,ib) = logSumExp([Pxy(jgb,ib);Pxy(jgx,ix)+Pxy(jxb,ib)-vnpxx(ix)]); % (8) in Wales2009, for b~=g
	endif
      else
	vn(ig)=vn(ig)+1; vn(ib)=vn(ib)+1;
	nn = size(Nxy,1);
	if or(vn(ig)>nn,vn(ib)>nn) Nxy = [Nxy; zeros(1,size(Nxy,2))]; Pxy = [Pxy; -99e9*ones(1,size(Pxy,2))]; endif
	Nxy(vn(ig),ig) = ib;
	Nxy(vn(ib),ib) = ig;
	jxg = find(Nxy(:,ig)==ix);
	jbg = vn(ig);
	jgb = vn(ib);
	if (vpxx(ix)<-0.5)
	  Pxy(jgb,ib) = logSumExp([Pxy(jgb,ib);Pxy(jgx,ix)+Pxy(jxb,ib)-logDiffExp(0,vpxx(ix))]); % (8) in Wales2009, for b~=g
	else
	  Pxy(jgb,ib) = logSumExp([Pxy(jgb,ib);Pxy(jgx,ix)+Pxy(jxb,ib)-vnpxx(ix)]); % (8) in Wales2009, for b~=g
	endif
      endif
    endfor
  endfor
% delete all connections to ix
  for jbx=1:nnow
    ib = vnow(jbx);
    jxb = find(Nxy(:,ib)==ix);
    Nxy(jxb,ib) = Nxy(vn(ib),ib);
    Nxy(vn(ib),ib) = 0;
    Pxy(jxb,ib) = Pxy(vn(ib),ib);
    Pxy(vn(ib),ib) = -99e9;
    vn(ib)--;
  endfor
endif
Nxy(:,ix) = [];
Pxy(:,ix) = [];
vpxx(ix) = [];
vnpxx(ix) = [];
vtaux(ix) = [];
vn(ix) = [];
n = size(vn,1);
% reindex
for ib = 1:n
  if (vn(ib)>0)
    for jxb = 1:vn(ib)
      if (Nxy(jxb,ib)>ix)
	Nxy(jxb,ib)--;
      endif
    endfor
  endif
endfor
endfunction