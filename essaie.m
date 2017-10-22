%% Pour tester la génération de loi uniforme sur {0,..,n-1} en utilisant la méthode générale
n=100;
n_sim=2000;

X=zeros(n_sim,1);

for i=1:n_sim
    X(i)=uniformGen(n);
end

moy=mean(X);
var=std(X)^2;


%% test de la fonction monteCarlo
clc
clear

g=@(x)log(x);
f=@(x) 1/pi;
[I_hat,flag,err_std,nb_sim]=monteCarlo(0.95,0.01,10000000,g,f);

%% test de comment enregistrer automatiquement les figures
x=0:0.01:pi;
figure(2);
plot(x)
chem='images';
chem=strcat(chem,'/droite0_pi');
print('-f2',chem,'-dpng')

figure(2);
plot(sin(x))
print('-f2','sin0_pi','-dpng')

%% test de comment enregistrer automatiquement les fichiers CSV
clc
clear

%delete test.csv
A=[23,24,67,4;87,13,999,7;65,67,546,13];
data=num2cell(A,2);
%c = {'abec' 'deddf' 'ghk' 'ggas18';23,24,67,4;87,13,999,7;65,6767,546,13};
c = cell(size(A,1)+1,size(A,2));
c(1,:)={'abec' 'deddf' 'ghk' 'ggas18'};
c(2:end,:)=num2cell(A);
fid = fopen('test.csv', 'w') ;
fprintf(fid, '%s,', c{1,1:end-1}) ;
fprintf(fid, '%s\n', c{1,end}) ;
fclose(fid) ;

dlmwrite('test.csv', c(2:end,:), '-append') ;