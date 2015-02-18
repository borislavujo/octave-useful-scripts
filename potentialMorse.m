function [vEne, Grad] = potentialMorse(X, nat, epsilon, range)
% calculation of potential and gradient for k replicas of a Morse cluster of nat atoms
% X is k*nat x 3 matrix of coordinates
% particle size is 1, default values for range and depth of the potential are 30 and 1, respectively
if (nargin<2) nat = 6; endif % number of atoms in each cluster
if (nargin<3) epsilon = 1; endif % epsilon
if (nargin<4) range = 30; endif % rho
n = size(X,1); % number of atoms altogether - each line of X corresponds to 1 atom
k = n/nat; % number of non-interacting Morse clusters in our sample

% first we need to make an array of distance matrices
XX = repmat(X(:,1),1,nat);
YY = repmat(X(:,2),1,nat);
ZZ = repmat(X(:,3),1,nat);
XXT = zeros(size(XX)); YYT=XXT; ZZT=XXT;
for i=1:k
  XXT((i-1)*nat+1:i*nat,:) = XX((i-1)*nat+1:i*nat,:)';
  YYT((i-1)*nat+1:i*nat,:) = YY((i-1)*nat+1:i*nat,:)';
  ZZT((i-1)*nat+1:i*nat,:) = ZZ((i-1)*nat+1:i*nat,:)';
endfor
DX = XX-XXT; DY = YY-YYT; DZ = ZZ-ZZT;
R = sqrt(DX.^2+DY.^2+DZ.^2); % array of dist matrices (size n x nat)
R(R==0)=9e99; % a dirty trick to avoid infinities at the diagonal

% calculation of energies for each cluster
T = exp(range*(ones(n,nat)-R)); % auxiliary variable
E = epsilon*T.*(T-2*ones(n,nat)); % array of atom-atom interaction energies
vE = sum(E,2); % each column is the interaction energy of particular atom with the rest in that cluster
vEne = zeros(k,1);
for i=1:k
  vEne(i) = sum(vE((i-1)*nat+1:i*nat))/2;
endfor

% calculation of the gradient
G = -2*epsilon*range*T.*(T-ones(n,nat)); % common term in all gradient components
GX = G.*DX./R; GY = G.*DY./R; GZ = G.*DZ./R;
Grad = [sum(GX,2), sum(GY,2), sum(GZ,2)];
