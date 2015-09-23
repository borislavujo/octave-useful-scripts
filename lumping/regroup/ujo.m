R = load('R').R;
% avoid grouping observables
va{1} = [1:5]'; va{2} = [401:5310]'; va{3} = [6:400]';
for i1=1:2
  for i2=i1+1:3
    R(va{i1},va{i2}) *= 1e-30;
    R(va{i2},va{i1}) *= 1e-30;
  endfor
endfor
vp = load('vp').vp;
[K,vi] = regroupfree(R,vp,1e-30,1210);
vj = unique(vi);
vi0 = vi; nFin = size(vj,1)
for i=1:size(vj,1)
  vi(vi0==vj(i))=i;
endfor
%vi = vi';
save K K
save vi vi