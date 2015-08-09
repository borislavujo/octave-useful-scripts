vTemps = [3:0.5:9]'; % vector of betas (inverse temperatures)
nSteps = 100000; % numer of steps for each temperature
dh = 0.01; % step size for MC
% load starting structure
nCells = 37; % 37 disctinct structures for Morse6(2D)
nTemp = size(vTemps,1); vInds=[1:nTemp]';
X = center(load('X')); [nAt,nDim]=size(X(:,1:2))
vEne = repmat(potential(X),nTemp,1);
X = repmat(X,nTemp,1);
Pops = zeros(nTemp,nCells);
Enes = zeros(nTemp,50); % also check energy distribution from -9 to -4
% initialise output arrays
acceptance = 0;
for iStep=1:nSteps
  showStep = iStep
  for iStr=1:nTemp
    X0 = X((iStr-1)*nAt+1:iStr*nAt,:);
    [X0,dene,acc] = mcmorse(X0,dh,vTemps(vInds(iStr)));
    iCell = wcell(X0);
    Pops(vInds(iStr),iCell)+=1;
    X((iStr-1)*nAt+1:iStr*nAt,:) = X0;
    vEne(iStr) += dene;
    rowInd = min(ceil((vEne(iStr)+9.1)/0.1),50);
    Enes(vInds(iStr),rowInd)+=1;
    acceptance += acc/(nSteps*nTemp);
  endfor
%  vstrs = randperm(nTemp)(1:2); s1=vstrs(1); s2=vstrs(2); % random pair
  s1 = ceil(rand(1)*(nTemp-1)); s2 = s1+1; % random neighbours
  expDEDb = exp((vEne(s1)-vEne(s2))*(vTemps(vInds(s1))-vTemps(vInds(s2))))
  if (expDEDb>rand(1))
    exchanging = [s1,s2]
    ujo = vInds(s1); vInds(s1) = vInds(s2); vInds(s2) = ujo;
%    ShowPop = [Pops(:,1:15),sum(Pops(:,16:36),2),Pops(:,37)]
  endif
  X0 = X((s1-1)*nAt+1:s1*nAt,:);
  iCell = wcell(X0);
  Pops(vInds(s1),iCell)+=1;
  rowInd = min(ceil((vEne(s1)+9.1)/0.1),50);
  Enes(vInds(s1),rowInd)+=1;
  X0 = X((s2-1)*nAt+1:s2*nAt,:);
  iCell = wcell(X0);
  Pops(vInds(s2),iCell)+=1;
  rowInd = min(ceil((vEne(s2)+9.1)/0.1),50);
  Enes(vInds(s2),rowInd)+=1;
endfor

acc = acceptance
save Pops Pops
save Enes Enes