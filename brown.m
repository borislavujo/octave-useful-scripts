function X = brown(X, beta, gamma, dt)
% propagate positions X (each line represents one particle) by time dt using Brownian dynamics
[n,k] = size(X);
[ePot,Xttm] = potential(X); % calculation of energy and energy gradient
Xt = -Xttm/gamma + sqrt(2/(beta*gamma*dt))*normrnd(0,1,n,k);
X = X + dt*Xt;