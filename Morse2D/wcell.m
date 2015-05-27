function cellInd = wcell(X)
% finds the cell index
X = center(X);
[nAt,nDim]=size(X); % should be 6,2
R = sqrt((repmat(X(:,1),1,nAt)-repmat(X(:,1)',nAt,1)).^2+(repmat(X(:,2),1,nAt)-repmat(X(:,2)',nAt,1)).^2);
bonLen = 1.1; % absolutely sufficient
CM = (R<bonLen)*1.0 - eye(nAt);
vBonds = sum(CM,2);
v = flipud(sort(vBonds))'; % sorted contact orders
nBonds = sum(vBonds)/2;

% 9-bond clusters
if (nBonds==9)
  if prod(v==[4,4,4,2,2,2])
    cellInd = 1;
  elseif prod(v==[5,3,3,3,2,2])
    cellInd = 2;
  elseif prod(v==[4,4,3,3,2,2])
    if chirTest34(X,CM)
      cellInd = 3;
    else
      cellInd = 4;
    endif
  else cellInd = 40;
  endif
% 8-bond clusters
elseif (nBonds==8)
  if (v(1)==5)
    if prod(v==[5,3,2,2,2,2])
      cellInd = 13;
    elseif prod(v==[5,3,3,2,2,1])
      cellInd = 12;
    else cellInd = 39;
    endif
  elseif (v(1)==4)
    if prod(v==[4,4,3,2,2,1])
      if chirTest78(X,CM)
	cellInd = 7;
      else
	cellInd = 8;
      endif
    elseif prod(v==[4,3,3,3,2,1])
      if chirTest56(X,CM)
	cellInd = 5;
      else
	cellInd = 6;
      endif
    elseif prod(v==[4,3,3,2,2,2])
      if not(test33(CM))
	cellInd = 14;
      elseif test323(CM)
	cellInd = 9;
      elseif chirTest1011(X,CM)
	cellInd = 10;
      else
	cellInd = 11;
      endif
    else cellInd = 39;
    endif
  elseif prod(v==[3,3,3,3,2,2])
    cellInd = 15;
  else cellInd = 39;
  endif
% 7-bond clusters
elseif (nBonds==7)
  if (v(6)==0)
    cellInd = 36;
  elseif (v(1)==5)
    if prod(v==[5,3,2,2,1,1])
      cellInd = 33;
    elseif prod(v==[5,2,2,2,2,1])
      cellInd = 34;
    else cellInd = 38;
    endif
  elseif (v(1)==4)
    if prod(v==[4,3,2,2,2,1])
      if test12(CM)
	cellInd = 21;
      elseif test13(CM)
	if chirTest1920(X,CM)
	  cellInd = 19;
	else 
	  cellInd = 20;
	endif
      else
	if chirTest2728(X,CM)
	  cellInd = 27;
	else cellInd = 28; endif
      endif
    elseif prod(v==[4,3,3,2,1,1])
      if test141(CM)
	cellInd = 30;
      elseif chirTest1718(X,CM)
	cellInd = 17;
      else
	cellInd = 18;
      endif
    elseif prod(v==[4,4,2,2,1,1])
      cellInd = 29;
    elseif prod(v==[4,2,2,2,2,2])
      cellInd = 32;
    else cellInd = 38;
    endif
  elseif (v(1)==3)
    if prod(v==[3,3,3,3,1,1])
      cellInd = 24;
    elseif prod(v==[3,3,3,2,2,1])
      if test12(CM)
	cellInd = 16;
      elseif test22(CM)
	cellInd = 26;
      elseif chirTest2223(X,CM)
	cellInd = 22;
      else
	cellInd = 23;
      endif
    elseif prod(v==[3,3,2,2,2,2])
      if test222(CM)
	cellInd = 31;
      elseif test223(CM)
	cellInd = 25;
      else cellInd = 35;
      endif
    else cellInd = 38;
    endif
  else cellInd = 38;
  endif
% 6- and fewer bonds -> state 37
else
  cellInd = 37;
endif

endfunction

%
%  AUXILIARY FUNCTIONS
%

%
% CHIRALITY TESTS
%
function yesorno = chirTest34(X,C)
v = sum(C,2);
vi2 = find(v==2);
i2 = vi2(1); % take one of the furthest disks
vi = find(C(:,i2)); % find its contacts
if (v(vi(1))==4) i4 = vi(1); else i4 = vi(2); endif % find the neighbour w 4 contacts
vc = cross([X(i4,:),0],[X(i2,:),0]);
if (vc(3)<0) yesorno = 1; else yesorno = 0; endif
endfunction
%
function yesorno = chirTest56(X,C)
v = sum(C,2);
i1 = find(v==1);
i4 = find(v==4);
vc = cross([X(i4,:),0],[X(i1,:),0]);
if (vc(3)<0) yesorno = 1; else yesorno = 0; endif
endfunction
%
function yesorno = chirTest78(X,C)
v = sum(C,2);
i1 = find(v==1);
i3 = find(v==3);
vc = cross([X(i3,:),0],[X(i1,:),0]);
if (vc(3)<0) yesorno = 1; else yesorno = 0; endif
endfunction
%
function yesorno = chirTest1011(X,C)
v = sum(C,2);
i4 = find(v==4);
vi2 = find(C(:,i4));
i21 = vi2(1);
vtest = C(:,i21);
vtest(i4) = 0;
itest = find(vtest);
if (v(itest)==3) i22 = vi2(2); else i22 = vi2(1); i21 = vi2(2); endif
vc = cross([X(i21,:),0],[X(i22,:),0]);
if (vc(3)<0) yesorno = 1; else yesorno = 0; endif
endfunction
%
function yesorno = chirTest1718(X,C)
v = sum(C,2);
vi1 = find(v==1);
i11 = vi1(1);
itest = find(C(:,i11));
if (v(itest)==3) i12 = vi1(2); else i12 = vi1(1); i11 = vi1(2); endif
vc = cross([X(i11,:),0],[X(i12,:),0]);
if (vc(3)<0) yesorno = 1; else yesorno = 0; endif
endfunction
%
function yesorno = chirTest1920(X,C)
v = sum(C,2);
vi2 = find(v==2);
if (C(vi2(1),vi2(2))==0 && C(vi2(1),vi2(3))==0) i2 = vi2(1);
elseif (C(vi2(2),vi2(3))==0 && C(vi2(2),vi2(1))==0) i2 = vi2(2);
else i2 = vi2(3); endif
i1 = find(v==1);
vc = cross([X(i2,:),0],[X(i1,:),0]);
if (vc(3)<0) yesorno = 1; else yesorno = 0; endif
endfunction
%
function yesorno = chirTest2223(X,C)
v = sum(C,2);
i1 = find(v==1);
vi3 = find(v==3);
if (C(vi3(1),vi3(2))==1 && C(vi3(1),vi3(3))==1) i3 = vi3(1); vi3(1) = [];
elseif (C(vi3(2),vi3(3))==1 && C(vi3(2),vi3(1))==1) i3 = vi3(2); vi3(2) = [];
else i3 = vi3(3); vi3(3) = []; endif
ic = vi3(1);
vtest = C(:,ic); vtest(i3) = 0; vitest = find(vtest);
if (v(vitest(1))~=2 || v(vitest(2))~=2) ic = vi3(2); endif
vc = cross([X(i3,:),0]-[X(ic,:),0],[X(i1,:),0]-[X(ic,:),0]);
if (vc(3)<0) yesorno = 1; else yesorno = 0; endif
endfunction
%
function yesorno = chirTest2728(X,C)
v = sum(C,2);
ic = find(v==4);
i1 = find(v==1);
vi2 = find(v==2);
for i=1:3
  vtest = find(C(:,vi2(i)));
  if (v(vtest(1))==3 && v(vtest(2))==4) i2 = vi2(i); endif
  if (v(vtest(1))==4 && v(vtest(2))==3) i2 = vi2(i); endif
endfor
vc = cross([X(i2,:),0]-[X(ic,:),0],[X(i1,:),0]-[X(ic,:),0]);
if (vc(3)<0) yesorno = 1; else yesorno = 0; endif
endfunction
%
%  GRAPH PATTERN SEARCHES
%
% test for 1-2 pattern in a contact matrix
function yesorno = test12(C)
vBonds = sum(C,2);
vi1 = find(vBonds==1);
vi2 = find(vBonds==2);
yesorno = 0;
for i=1:size(vi2,1)
  if (C(vi1(1),vi2(i))==1)
    yesorno = 1;
  endif
endfor
endfunction
%
% test for 1-3 pattern in a contact matrix
function yesorno = test13(C)
vBonds = sum(C,2);
vi1 = find(vBonds==1);
vi3 = find(vBonds==3);
yesorno = 0;
for i=1:size(vi3,1)
  if (C(vi1(1),vi3(i))==1)
    yesorno = 1;
  endif
endfor
endfunction
%
% test for 1-4-1 pattern in a contact matrix
% in a matrix with contacts 433211
function yesorno = test141(C)
vBonds = sum(C,2);
vi1 = find(vBonds==1);
vi4 = find(vBonds==4);
yesorno = 0;
if (C(vi1(1),vi4(1))==1 && C(vi1(2),vi4(1))==1)
  yesorno = 1;
endif
endfunction 
%
% test for 2-2 pattern in a contact matrix
% with contacts 333221
function yesorno = test22(C)
vBonds = sum(C,2);
vi2 = find(vBonds==2);
yesorno = 0;
if (C(vi2(1),vi2(2))==1)
  yesorno = 1;
endif
endfunction
%
% test for 2-2-2 pattern in a contact matrix
% for a system with 332222 contacts
function yesorno = test222(C)
vBonds = sum(C,2);
vi2 = find(vBonds==2);
if (C(vi2(1),vi2(2))==1 && C(vi2(1),vi2(3))==1)
  yesorno = 1;
elseif (C(vi2(2),vi2(1))==1 && C(vi2(2),vi2(3))==1)
  yesorno = 1;
elseif (C(vi2(3),vi2(1))==1 && C(vi2(3),vi2(2))==1)
  yesorno = 1;
else
  yesorno = 0;
endif
endfunction
%
% test for 2-2-3 pattern in a contact matrix
% for a system with 332222 contacts
function yesorno = test223(C)
vBonds = sum(C,2);
vi3 = find(vBonds==3);
vi2 = find(vBonds==2);
yesorno = 0;
for i21=1:3
  for i22=i21+1:4
    for i3=1:2
      t22 = C(vi2(i21),vi2(i22));
      t321 = C(vi2(i21),vi3(i3));
      t322 = C(vi2(i22),vi3(i3));
      if (t22*t321*t322==1)
	yesorno = 1;
      endif
    endfor
  endfor
endfor
endfunction
%
% test for 3-3 pattern in a contact matrix
% for a system with 433222 contacts
function yesorno = test33(C)
vBonds = sum(C,2);
vi3 = find(vBonds==3);
if (C(vi3(1),vi3(2))==1)
  yesorno = 1;
else
  yesorno = 0;
endif
endfunction
%
% test for 3-2-3 pattern in a contact matrix
% for a system with 433222 contacts
function yesorno = test323(C)
vBonds = sum(C,2);
vi3 = find(vBonds==3);
vi2 = find(vBonds==2);
if (C(vi3(1),vi2(1))==1 && C(vi3(2),vi2(1))==1)
  whatswrong = 1
  yesorno = 1;
elseif (C(vi3(1),vi2(2))==1 && C(vi3(2),vi2(2))==1)
  yesorno = 1;
elseif (C(vi3(1),vi2(3))==1 && C(vi3(2),vi2(3))==1)
  yesorno = 1;
else
  yesorno = 0;
endif
endfunction
