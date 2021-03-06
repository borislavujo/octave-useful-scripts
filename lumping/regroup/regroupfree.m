function [R,vi] = regroupfree(R,vp,thresh,nDesiredStates)
% grouping minima stops when either a rate constant threshold or a number of states is reached
  if (nargin<4) nDesiredStates = 0; endif 
  n = size(R,1); nStates = n; % the original number of states
  vi = [1:n]'; % index vector, always of size n
  R = R - diag(diag(R)); % delete all diagonal (negative) terms
  [vmax,vrowi] = max(min(R,R')); [kmax,icol] = max(vmax); % find the largest smaller rate constant
  while and(kmax>thresh,nStates>nDesiredStates)
    [R,vp,vi] = groupStates(R,vi,vp,vrowi(icol),icol); % both vp and R are updated
    [vmax,vrowi] = max(min(R,R')); [kmax,icol] = max(vmax); % find the largest smaller rate constant
    know = full(kmax) % print progress
    nStates = size(R,1) % print progress
  endwhile
  R = R - diag(sum(R)); % regrouped rate matrix -> add diagonal terms
endfunction

function [R,vp,vi] = groupStates(R,vi,vp,is1,is2)
% rate matrix reduction by 1 state
  R(is1,is2) = 0; R(is2,is1) = 0; % delete rates between s1 and s2, so they do not appear among neighbours
  vnei = unique([find(R(is1,:)),find(R(is2,:))]); % find all the neighbours of s1 and s2
  p1 = vp(is1); p2 = vp(is2); % populations of s1 and s2
  for j=vnei
    R(is1,j) = (p1*R(is1,j)+p2*R(is2,j))/(p1+p2); % rate from a grouped state to a neighbour
    R(j,is1) = R(j,is1) + R(j,is2); % rate from a neighbour to a grouped state
  endfor
  vp(is1) = p1+p2; % grouped population
  vp(is2) = []; R(is2,:) = []; R(:,is2) = []; % delete s2 from R and vp
  vi([find(vi==is1);find(vi==is2)]) = is1; % indices of all minima grouped into s1
  nStates = size(R,1);
  if (is2<=nStates)
    vi(find(vi>=is2))--; % deleting s2, R became smaller -> decrement all indices>=is2
  endif
endfunction
