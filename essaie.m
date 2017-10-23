%% test de la fonction monteCarlo
clc
clear

g=@(x)log(x);
f=@(x) 1/pi;
[I_hat,flag,err_std,nb_sim]=monteCarlo0(0.95,0.01,10000000,g,f);

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

%% monte carlo pour l'exercice 1 call

clc
clear
close all

N=100000;
Nsim=500;
pas=10;
I=zeros(1,Nsim);
Err=I;
K=1;
bet=1;
C=(exp(bet^2/2)*normcdf(bet-log(K)/bet,0,1)-K*normcdf(-log(K)/bet,0,1))*ones(1,Nsim);
for n=1:Nsim
    
[I_hat,err_std]=monteCarloCall(0.95,n*pas);
I(n)=I_hat;
Err(n)=err_std;
end

% affichage des resultats
figure();
hold on
figEx=plot(C);
figI=plot(I);
figHaut=plot(I+Err);
figBas=plot(I-Err);
title('Monte carlo exo1 Call')
legend([figEx, figI, figHaut, figBas],'exact', 'estimation','borne haute','borne basse');

%% monte carlo pour l'exercice 1 put

clc
clear
close all

N=100000;
Nsim=500;
pas=100;
I=zeros(1,Nsim);
Err=I;
K=1;
bet=1;
C=(-exp(bet^2/2)*normcdf(-bet+log(K)/bet,0,1)+K*normcdf(log(K)/bet,0,1))*ones(1,Nsim);
for n=1:Nsim
    
[I_hat,err_std]=monteCarloPut(0.95,n*10);
I(n)=I_hat;
Err(n)=err_std;
end

% affichage des resultats
figure();
hold on
figEx=plot(C);
figI=plot(I);
figHaut=plot(I+Err);
figBas=plot(I-Err);
title('Monte carlo exo1 Put')
legend([figEx, figI, figHaut, figBas],'exact', 'estimation','borne haute','borne basse');


