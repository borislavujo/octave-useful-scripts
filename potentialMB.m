function [vEne, Grad] = potentialMB(X)
% energy and gradient for the Muller-Brown potential (m paricles in 2D)
A = [-200, -100, -170, 15]; % parameters of the MB potential
a = [-1, -1, -6.5, 0.7];
b = [0, 0, 11, 0.6];
c = [-10, -10, -6.5, 0.7];
x0 = [1, 0, -0.5, -1];
y0 = [0, 0.5, 1.5, 1.0];
m = size(X,1); % number of (non-interacting) particles - the calculation is vectorised
vEne = zeros(m,1); % array initiation
Grad = zeros(m,2);
% loop over all 4 contributions
for i=1:4
  XX = X(:,1)-x0(i); % auxiliary variable
  YY = X(:,2)-y0(i);
  vEne = vEne + A(i)*exp(a(i)*XX.*XX + b(i)*XX.*YY + c(i)*YY.*YY); % energy
  Grad(:,1) = Grad(:,1) + A(i)*(2*a(i)*XX+b(i)*YY).*exp(a(i)*XX.*XX + b(i)*XX.*YY + c(i)*YY.*YY); 
  Grad(:,2) = Grad(:,2) + A(i)*(b(i)*XX+2*c(i)*YY).*exp(a(i)*XX.*XX + b(i)*XX.*YY + c(i)*YY.*YY);
endfor
