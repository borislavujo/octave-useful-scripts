function [ePot, Grad] = potentialLJ(X);
% calculates potential energy and gradient of a Lennard-Jones cluster
% in 1D, 2D or 3D
[nAt,nDim] = size(X);
vx = X(:,1);
if (nDim>1) vy=X(:,2); else vy=zeros(nAt,1); endif
if (nDim>2) vz=X(:,3); else vz=zeros(nAt,1); endif

Dx = repmat(vx,1,nAt) - repmat(vx',nAt,1); Dy = repmat(vy,1,nAt) - repmat(vy',nAt,1); Dz = repmat(vz,1,nAt) - repmat(vz',nAt,1);
R12s = Dx.*Dx + Dy.*Dy + Dz.*Dz;
R12 = sqrt(R12s);
Orientx = Dx ./ (R12+eye(nAt)); Orienty = Dy ./ (R12+eye(nAt)); Orientz = Dz ./ (R12+eye(nAt));
R6term = R12s.^3;
R12recip = ones(nAt)./(R12 + eye(nAt)) - eye(nAt);
R6recip = ones(nAt)./(R6term + eye(nAt)) - eye(nAt);
Forces = -24*(2*R6recip - ones(nAt)).*R6recip.*R12recip;
Forcex = Orientx .* Forces; Forcey = Orienty .* Forces; Forcez = Orientz .* Forces;
vfx = Forcex * ones(nAt,1); vfy = Forcey * ones(nAt,1); vfz = Forcez * ones(nAt,1);
if (nDim==1) Grad = vfx;
elseif (nDim==2) Grad = [vfx,vfy];
else Grad = [vfx, vfy, vfz]; endif

PotE = 4*R6recip.*(R6recip-ones(nAt));
ePot = 0.5*(4*sum(sum(R6recip.*(R6recip-ones(nAt)))));