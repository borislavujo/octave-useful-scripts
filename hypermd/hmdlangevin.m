function [X, Xt, dtScaled] = hmdlangevin(X, Xt, beta, gamma, dt, vm, epot0, alpha)
% propagate positions X and momenta Xt (each line represents one particle) by time dt
% one replica with n particles in k dimensions
% implementation of OVRVO integrator adapted from Sivak et al., 2013, arXiv 1301.3800
% biasing potential according to Hamelberg et al., 2004, JCP 120, 11919
[n,k] = size(X);
[n2,k2] = size(Xt);
% sanity check of input arrays
if (n~=n2) error("X and Xt must have same sizes (number of particles)"); endif
if (k~=k2) error("X and Xt must have same sizes (dimension of space)"); endif
if (nargin<6) vm = ones(n,1); endif % assign default masses
if (size(vm,1)~=n) error("number of masses not equal to number of particles"); endif

% a, b ... auxiliary variables
va = exp(-gamma*dt*ones(n,1)./vm);
if (gamma == 0) vb = ones(n,1); else vb = sqrt(tanh(0.5*gamma*dt*ones(n,1)./vm)./(0.5*gamma*dt*ones(n,1)./vm)); endif % limit b (gamma=0) = 1 

% OVRVO integrator, equations 7a-g in arXiv 1301.3800
Xt = repmat(sqrt(va),1,k).*Xt + repmat(sqrt(ones(n,1)-va),1,k).*normrnd(0,1,n,k)./sqrt(beta*repmat(vm,1,k)); % 7a 
[epot,Xttm] = potential(X); % calculation of energy and energy gradient
epotb = (epot0-epot).^2./(alpha+epot0-epot);
dtScaled = exp(beta*epotb)*dt/2;
Xttm = Xttm * (1 + (epot0-epot)^2/(alpha+epot0-epot)^2 + 2*(epot0-epot)/(alpha+epot0-epot);
Xt = Xt - ((dt/2)*repmat(vb,1,k).*Xttm./repmat(vm,1,k)); % 7b (acceleration = - energy gradient / mass)
X = X + (dt/2)*repmat(vb,1,k).*Xt; % 7c
[epot,Xttm] = potential(X); % 7d evaluate gradient at midpoint
epotb = (epot0-epot).^2./(alpha+epot0-epot);
dtScaled = dtScaled + exp(beta*epotb)*dt/2;
Xttm = Xttm * (1 + (epot0-epot)^2/(alpha+epot0-epot)^2 + 2*(epot0-epot)/(alpha+epot0-epot);
X = X + ((dt/2)*repmat(vb,1,k).*Xt); % 7e
Xt = Xt - ((dt/2)*repmat(vb,1,k).*Xttm./repmat(vm,1,k)); % 7f
Xt = repmat(sqrt(va),1,k).*Xt + sqrt(repmat(ones(n,1)-va,1,k)).*normrnd(0,1,n,k)./sqrt(beta*repmat(vm,1,k)); % 7g
