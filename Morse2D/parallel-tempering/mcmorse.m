function [X, dene, accepted] = mcmorse(X, dh, beta, epsilon, range, rmax)
% updates positions of a morse cluster according to Metropolis criterion
% X is n x 3 matrix of coordinates
if (nargin<2) dh = 0.01; endif % epsilon
if (nargin<3) beta = 6; endif % rho
if (nargin<4) epsilon = 1; endif % epsilon
if (nargin<5) range = 30; endif % rho
[nAt,nDim] = size(X);
if (nargin<6) rmax = 3; endif % max distance from the centre
% randomly select a particle
ipart = ceil(nAt*rand(1));
vx = X(ipart,:);
Xrest = X; Xrest(ipart,:) = []; % remaining nAt-1 atoms
% energy now
Dxyz0 = repmat(vx,nAt-1,1)-Xrest;
vr0 = sqrt(sum(Dxyz0.^2,2));
vt0 = exp(range*(ones(nAt-1,1)-vr0));
ene0 = sum(epsilon*vt0.*(vt0-2*ones(nAt-1,1)));
% move
% check for geometric criteria first
mrok = 0;
while (mrok==0)
  vxn = vx + normrnd(0,dh,1,nDim);
  XX = [Xrest;vxn];
  XX = center(XX); % center
  mr = max(sum(XX.*XX,2)); % maximum square distance
  if (mr<rmax^2) mrok=1; endif
endwhile
% calculate the energy change
Dxyz = repmat(vxn,nAt-1,1)-Xrest;
vr = sqrt(sum(Dxyz.^2,2));
vt = exp(range*(ones(nAt-1,1)-vr));
ene = sum(epsilon*vt.*(vt-2*ones(nAt-1,1)));
% accepted?
accepted = 0;
dene = 0;
if (rand(1)<exp(beta*(ene0-ene)));
  accepted = 1;
  dene = ene-ene0;
  X(ipart,:) = vxn;
endif
% center
X = center(X);