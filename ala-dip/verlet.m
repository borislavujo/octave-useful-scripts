function [X, Xt, Grad] = verlet(X, Xt, dt, vm, Grad)
% Velocity Verlet algorithm, propagates positions (X) and velocities (Xt) by time dt
[nAt,nDim] = size(X); [nat,ndim] = size(Xt);
if (nat~=nAt | ndim~=nDim) error('dimensions of X and Xt disagree'); endif % check dimensions
if (nargin<4) vm = ones(nAt,1); endif % set default masses
% if force in previous point is not given, calculate gradient
if (nargin<5) [ene, Grad] = potentialAMB(X); endif

Xt = Xt - 0.5*dt*Grad./repmat(vm,1,nDim); % 1. update velocities
X = X + dt*Xt; % 2. update positions
[ene, Grad] = potentialAMB(X); % 3. calculate gradient
Xt = Xt - 0.5*dt*Grad./repmat(vm,1,nDim); % 4. update velocities
eneKin = sum(sum(Xt.^2,2).*vm)/2
enePot = ene
eneTot = eneKin + enePot
