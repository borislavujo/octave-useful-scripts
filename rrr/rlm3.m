function [Kl, vpl] = rlm3(m1,m2)
% returns 3x3 log matrix, diagonal is log(-K(i,i))
if (nargin<2) m2 = 200; endif % magnitude of the equil pop distribution
if (nargin<1) m1 = 200; endif % magnitude of the rxn time distribution
vpl = m2*rand(3,1); vpl = vpl - logSumExp(vpl);
vtaul = m1*rand(3,1);
vtaul = [min(vtaul); vtaul(find(vtaul>min(vtaul)))];
Kl = zeros(3);
Kl(1,2) = vpl(1) - logSumExp([vpl(1);vpl(2)]) - vtaul(1);
Kl(1,3) = vpl(1) - logSumExp([vpl(1);vpl(3)]) - vtaul(2);
Kl(2,3) = vpl(2) - logSumExp([vpl(2);vpl(3)]) - vtaul(3);
Kl(2,1) = Kl(1,2) + vpl(2) - vpl(1);
Kl(3,1) = Kl(1,3) + vpl(3) - vpl(1);
Kl(3,2) = Kl(2,3) + vpl(3) - vpl(2);
Kl(1,1) = logSumExp(Kl(2:3,1));
Kl(2,2) = logSumExp(Kl([1;3],2));
Kl(3,3) = logSumExp(Kl(1:2,3));
endfunction