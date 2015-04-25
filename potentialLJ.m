function [ePot, Grad] = potentialLJ(X);
% calculates potential energy and gradient of a Lennard-Jones cluster
% in 1D, 2D or 3D
[nAt,nDim] = size(vx);
vx = X(:,1); 
if (nDim>1) vy=X(:,2); else vy=zeros(nAt,1); endif
if (nDim>2) vz=X(:,3); else vz=zeros(nAt,1); endif

Dx = repmat(vx,1,n) - repmat(vx',n,1); Dy = repmat(vy,1,n) - repmat(vy',n,1); Dz = repmat(vz,1,n) - repmat(vz',n,1);
R12s = Dx.*Dx + Dy.*Dy + Dz.*Dz;
R12 = sqrt(R12s);
Orientx = Dx ./ (R12+eye(n)); Orienty = Dy ./ (R12+eye(n)); Orientz = Dz ./ (R12+eye(n));
R6term = R12s.^3;
R12recip = ones(n)./(R12 + eye(n)) - eye(n);
R6recip = ones(n)./(R6term + eye(n)) - eye(n);
Forces = -24*(2*R6recip - ones(n)).*R6recip.*R12recip;
Forcex = Orientx .* Forces; Forcey = Orienty .* Forces; Forcez = Orientz .* Forces;
vfx = - Forcex * ones(n,1); vfy = - Forcey * ones(n,1); vfz = - Forcez * ones(n,1);
if (nDim==1) Grad = vfx;
elseif (nDim==2) Grad = [vfx,vfy];
else Grad = [vfx, vfy, vfz]; endif

PotE = 4*R6recip.*(R6recip-ones(n));
ePot = 0.5*(4*sum(sum(R6recip.*(R6recip-ones(n)))));