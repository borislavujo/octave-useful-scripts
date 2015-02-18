function RotMat = ranRotMat(nDim)
% returns a random rotation matrix of size nDim x nDim
% nDim can be 2 or 3, default dimensionality is 3
% for any vector vx in 3D, RotMat*vx is randomly distributed on a sphere with radius sqrt(vx'*vx)
if (nargin<1) nDim = 3; endif
% random angle
  phi = 2*pi*rand(1);
  cp = cos(phi); sp = sin(phi);
  RR = [ cp, -sp ; ...
         sp,  cp ];
if (nDim==2)
  RotMat = RR;
elseif (nDim==3)
% method suggested by Arvo
  R = [RR, [0; 0] ; ...
        0,  0,  1 ];
  psi = 2*pi*rand(1);
  x = rand(1);
  vv = [cos(psi)*sqrt(x);sin(psi)*sqrt(x);sqrt(1-x)];
  H = eye(3) - 2*vv*(vv');
% result
  RotMat = H*R;
else
  error("argument of this function can be either 2 or 3");
endif