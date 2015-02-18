function [RotMat, dist] = qtrn3D(X,Y)
% finds the optimum rotation matrix RotMat (3 x 3), so that abs((X*RotMat - Y)^2) is minimised
% X ... n x 3, Y ... n x 3, RotMat ... 3 x 3, dist ... scalar
% adapted from Kearsley 1989, ActaCryst A45, 208
[m,n] = size(X);
% checking the input
[my,ny] = size(Y);
if (n~=3 || ny~=3)  printf("method is only for vectors in 2 dimensions");  return; 
elseif (my~=m) orubtf("mismatched number of atoms");  return; endif

% S is a 3x3 matrix
S = Y'*X;
% expression for K (size 4 x 4) is quite complicated
A = eye(2); B = fliplr(A); C = [1,0;0,-1]; D = fliplr(C); O = zeros(2); % auxiliary 2 x 2 matrices
K = S(1,1)*[A,O;O,-A] + S(1,2)*[O,B;B,O]  + S(1,3)*[O,-C;-C,O] + ...
    S(2,1)*[O,-D;D,O] + S(2,2)*[C,O;O,C]  + S(2,3)*[B,O;O,B]   + ...
    S(3,1)*[O,A;A,O]  + S(3,2)*[-B,O;O,B] + S(3,3)*[C,O;O,-C]  ;
% calculation of quaternions
[V,lam] = eig(K);
[smallestLam,indSmall] = min(diag(lam));
vq = V(:,indSmall);
% quaternions to rotation matrix
F = A*vq(1)+D*vq(4); G = C*vq(3)-B*vq(2); H = A*vq(3)-D*vq(2); % auxiliary 2 x 2 matrices
RotMatExt = [F,-G;G,F]*[F,-H';H,F'];
RotMat = -RotMatExt(1:3,1:3)';
dist = sum(sum((X*RotMat-Y).^2))/m;