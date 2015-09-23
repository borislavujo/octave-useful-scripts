R = load('R').R;
vp = load('vp').vp;
[K,vi] = regroupfree(R,vp,1e-5);
vj = unique(vi); vi0 = vi; nFin = size(vj,1)
for i=1:size(vj,1)
  vi(vi0==vj(i))=i;
endfor
vi = vi'
save K K
save vi vi