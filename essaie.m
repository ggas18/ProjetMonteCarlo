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

% le quantile d'ordre (alpha+1)/2 de la loi normale
% centrée reduite.
alpha=.95;
Z=norminv((alpha+1)/2,0,1); 
    
Nsim=500;
pas=10;
I=zeros(1,Nsim);
Err=I;
K=1;
beta=1;
C=(exp(beta^2/2)*normcdf(beta-log(K)/beta,0,1)-K*normcdf(-log(K)/beta,0,1))*ones(1,Nsim);
for n=1:Nsim
    
[I_hat,err_std]=monteCarloCall(n*pas);
I(n)=I_hat;
Err(n)=err_std;
end

% affichage des resultats
figure();
hold on
figEx=plot(C);
figI=plot(I);
figHaut=plot(I+Z*Err);
figBas=plot(I-Z*Err);
title('Monte carlo exo1 Call')
legend([figEx, figI, figHaut, figBas],'exact', 'estimation','borne haute','borne basse');

%% monte carlo pour l'exercice 1 put

clc
clear
close all

% le quantile d'ordre (alpha+1)/2 de la loi normale
% centrée reduite.
alpha=.95;
Z=norminv((alpha+1)/2,0,1); 

Nsim=500;
pas=10;
I=zeros(1,Nsim);
Err=I;
K=1;
beta=1;
P=(-exp(beta^2/2)*normcdf(-beta+log(K)/beta,0,1)+K*normcdf(log(K)/beta,0,1))*ones(1,Nsim);
for n=1:Nsim
    
[I_hat,err_std]=monteCarloPut(n*10);
I(n)=I_hat;
Err(n)=err_std;
end

% affichage des resultats
figure();
hold on
figEx=plot(P);
figI=plot(I);
figHaut=plot(I+Z*Err);
figBas=plot(I-Z*Err);
title('Monte carlo exo1 Put')
legend([figEx, figI, figHaut, figBas],'exact', 'estimation','borne haute','borne basse');

%% test distance de Levenshtein
clear all

s = char('tests');
t = char('twist');

% Edit Matrix
m=length(s);
n=length(t);
mat=zeros(m+1,n+1);
for i=1:1:m
    mat(i+1,1)=i;
end
for j=1:1:n
    mat(1,j+1)=j;
end
for i=1:m
    for j=1:n
        if (s(i) == t(j))
            mat(i+1,j+1)=mat(i,j);
        else
            mat(i+1,j+1)=1+min(min(mat(i+1,j),mat(i,j+1)),mat(i,j));
        end
    end
end

% Edit Sequence
s = char('tests');
t = char('twist');
i = m+1;
j = n+1;
display([s ' --> ' t])
while(i ~= 1 && j ~= 1)
    temp = min(min(mat(i-1,j-1), mat(i,j-1)), mat(i-1,j));
    if(mat(i-1,j) == temp)
        i = i - 1;
        t = [t(1:j-1) s(i) t(j:end)];
        disp(strcat(['insertion:' s ' --> ' t]))
    elseif(mat(i-1,j-1) == temp)
        if(mat(i-1,j-1) == mat(i,j))
            i = i - 1;
            j = j - 1;
            disp(strcat(['unchanged:' s ' --> ' t]))
        else
            i = i - 1;
            j = j - 1;
            t(j) = s(i);
            disp(strcat(['substition:' s ' --> ' t]))
        end
    elseif(mat(i,j-1) == temp)
        j = j - 1;
        t(j) = [];
        disp(strcat(['deletion:' s ' --> ' t]))
    end
end

%% optimisation des codes de monte carlo pour le call
clc
clear
close all

N_t=500;
pas=10;

[I_hat,IC_h,IC_b,C_exact]=mcCall(pas,N_t);
% affichage des resultats
figure();
hold on
figEx=plot(C_exact);
figI=plot(I_hat);
figHaut=plot(IC_h);
figBas=plot(IC_b);
title('Monte carlo exo1 Call')
legend([figEx, figI, figHaut, figBas],'exact', 'estimation','borne haute','borne basse');

%% optimisation des codes de monte carlo pour le call
clc
clear
close all

N_t=2000;
pas=1000;

[I_hat,IC_h,IC_b,P_exact]=mcPut(pas,N_t);
% affichage des resultats
figure();
hold on
figEx=plot(P_exact);
figI=plot(I_hat);
figHaut=plot(IC_h);
figBas=plot(IC_b);
title('Monte carlo exo1 Put')
legend([figEx, figI, figHaut, figBas],'exact', 'estimation','borne haute','borne basse');
