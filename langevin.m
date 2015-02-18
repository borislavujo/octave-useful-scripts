function [X, Xt] = langevin(X, Xt, beta, gamma, dt, vm)
% propagate positions X and momenta Xt (each line represents one particle) by time dt
% implementation of OVRVO integrator adapted from Sivak et al., 2013, arXiv 1301.3800
[n,k] = size(X);
[n2,k2] = size(Xt);
if (n~=n2) error("X and Xt must have same sizes (number of particles)"); endif
if (k~=k2) error("X and Xt must have same sizes (dimension of space)"); endif
if (nargin<6) vm = ones(n,1); endif % assign default masses
if (size(vm,1)~=n) error("number of masses not equal to number of particles"); endif

a = exp(-gamma*dt);
if (gamma == 0) 
  b = 1;
else
  ujo = gamma*dt/2;
  b = sqrt(tanh(ujo)/ujo);
endif
% OVRVO integrator
% 7a
Xt = sqrt(a)*Xt + sqrt(1-a)*normrnd(0,1,size(X,1),2)./sqrt(beta*repmat(vm,1,k));
% 7b
[vEne,Xttm] = potential(X); % calculation of energy and energy gradient
Xt = Xt - (dt/2)*b*Xttm./repmat(vm,1,k); % acceleration = - energy gradient / mass
% 7c-e
X = X + (dt/2)*b*Xt;
[vEne,Xttm] = potential(X); % evaluate gradient at midpoint
X = X + (dt/2)*b*Xt;
% 7f
Xt = Xt - (dt/2)*b*Xttm./repmat(vm,1,k);
% 7g
Xt = sqrt(a)*Xt + sqrt(1-a)*normrnd(0,1,size(X,1),2)./sqrt(beta*repmat(vm,1,k));
