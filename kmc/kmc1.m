function [indNew, t] = kmc1(ind, vN, IN, P, K, vK)
vI = IN(1:vN(ind),ind);
vP = P(1:vN(ind),ind);
randNum = rand(1);
vOK = vI(vP>randNum);
indNew = min(vOK);
t = -log(rand(1))/vK(ind);