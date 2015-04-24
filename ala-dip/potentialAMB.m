function [ene, Grad] = potentialAMB(X)
subor = fopen('start', 'w');
for i=1:size(X,1)
	fprintf(subor,"%20.9f%20.9f%20.9f\n",X(i,1),X(i,2),X(i,3));
endfor
system("cp start backup_pociatek");
fclose(subor)
%save -ascii start X ;
system("rm GRADIENTS ujo.out");
isok = system("./A9OPTIM > ujo.out");
Grad = load('GRADIENTS');
[statut,energy] = system("grep kcal ujo.out | awk '{printf $6}'");
ene = str2double(energy);
system("mv start pociatek");
