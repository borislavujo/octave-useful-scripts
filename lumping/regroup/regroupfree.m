function [R,vi] = regroupfree(R0,vp0,thresh)
n = size(R0,1);
R = sparse(R0);
vi = [1:n]';
vp = vp0;
R = R - diag(diag(R));
[vmax,vrowi] = max(min(R,R'));
[kmax,icol] = max(vmax);
while (kmax>thresh)
  [R,vp,vi] = groupStates(R,vi,vp,vrowi(icol),icol); % both vp and R are updated - info for group is the same for all elements in that group (so no change of matrix size / re-indexing) is needed
  [vmax,vrowi] = max(min(R,R'));
  [kmax,icol] = max(vmax);
endwhile
R = full(R); R = R - diag(sum(R));
endfunction

function [R,vp,vi] = groupStates(R,vi,vp,irow,icol)
is1 = vi(irow); is2 = vi(icol);
newind = min(is1,is2);
vi1 = find(vi==is1); vi2 = find(vi==is2);
R(vi1,vi2) = 0; R(vi2,vi1) = 0;
vnei = unique(vi(find(any(R([vi1;vi2],:)))));
p1 = vp(irow); p2 = vp(icol);
for j=vnei'
  vj = find(vi==j);
  R(vj,[vi1;vi2]) = (p1*R(vj(1),vi1(1))+p2*R(vj(1),vi2(1)))/(p1+p2);
  R([vi1;vi2],vj) = R(vi1(1),vj(1)) + R(vi2(1),vj(1));
endfor
vi([vi1;vi2]) = newind;
vp([vi1;vi2]) = p1+p2;
endfunction

