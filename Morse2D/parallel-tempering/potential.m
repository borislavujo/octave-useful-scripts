function [energy, Grad] = potential(X, epsilon, range)
% calculates potential energy
% X is n x 2 matrix of coordinates
if (nargin<2) epsilon = 1; endif % epsilon
if (nargin<3) range = 30; endif % rho
% particle size is 1
[n,k] = size(X); % n ... no of atoms; k .. dimension (2 or 3)
DX = repmat(X(:,1),1,n)-repmat(X(:,1)',n,1);
DY = repmat(X(:,2),1,n)-repmat(X(:,2)',n,1);
if (k>2) DZ = repmat(X(:,3),1,n)-repmat(X(:,3)',n,1); else DZ = zeros(size(DY)); endif
R = sqrt(DX.^2+DY.^2+DZ.^2);
R(R==0)=9e99;
T = exp(range*(ones(n)-R));
E = epsilon*T.*(T-2*ones(n));
E = E - tril(E);
energy = sum(sum(E));
G = -2*epsilon*range*T.*(T-ones(n));
GX = G.*DX./R;
GY = G.*DY./R;
GZ = G.*DZ./R;
if (k>2)
  Grad = [sum(GX')', sum(GY')', sum(GZ')'];
else
  Grad = [sum(GX')', sum(GY')'];
endif