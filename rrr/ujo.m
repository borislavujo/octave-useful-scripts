Kl = load('Kl').Kl
vpl = load('vpl')
vil = load('vlump');
vil(find(vil==0))=2
[vi,FEl,SSl,TSTl,NFl] = rrr(Kl,vpl,0,vil) % group into a 2-state model with arbitrary observables
