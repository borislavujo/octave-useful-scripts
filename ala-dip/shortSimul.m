%beta = 1.6778;
beta = 10;
%gamma = 55; % in reduced units, corresponds to ~1000 ps^{-1}
gamma = 0; % in reduced units, corresponds to ~100 ps^{-1}
dt0 = 0.020457; % 1 fs
nSteps = 1000;
dt = dt0/5;

% read the initial coordinates
X = load('pociatek');
[n,m] = size(X); % n .. no of atoms, m .. no of dimensions (3)
Xt = zeros(n,m);
vm = load('masses');
for i=1:n
    Xt(i,:) = normrnd(0,sqrt(1/(beta*vm(i))),1,m);
endfor

[ene,Grad] = potentialAMB(X);
Angles = zeros(nSteps,2);
for i=1:nSteps
%    dt = dtmax*min(i,linInc)/linInc; % linearly increase for linInc steps
%    dt = 1e-3;
    Phi = X([5,7,9,15]',:); phi = calcDihedral(Phi)
    Psi = X([7,9,15,17]',:); psi = calcDihedral(Psi)
    [X, Xt, Grad] = verlet(X, Xt, dt, vm, Grad);
%    langevin(X,Xt,beta,gamma,dt,vm);
    Angles(i,:) = [phi,psi];
endfor

save Angles Angles