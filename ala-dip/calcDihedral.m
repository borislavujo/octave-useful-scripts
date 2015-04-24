function dih = calcDihedral(X)
% X is a 4x3 matrix of positions (one atom per line)
% 1. unit vectors of all distances
VD = diff(X);
VU = VD ./ repmat(sqrt(sum(VD.^2,2)),1,3);
v1 = VU(1,:)'; v2 = VU(2,:)'; v3 = VU(3,:)';
va = v1 - (v1'*v2)*v2; 
vb = v3 - (v3'*v2)*v2;
va = va / sqrt(sum(va.^2));
vb = vb / sqrt(sum(vb.^2));
dih = pi-acos(va'*vb);
dih = dih * 360/(2*pi); % in degrees
if (v2'*cross(va,vb)>0) dih = -dih; endif
