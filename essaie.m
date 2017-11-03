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

%% monte carlo pour l'exercice 1

clc
close all

% le quantile d'ordre (alpha+1)/2 de la loi normale
% centree reduite.
alpha=.95;
Z=norminv((alpha+1)/2,0,1); 

% le nombre d'estimations des options.
Nsim=100;
% les pas entre la taille des echantillons
pas=2000;

% on initialise le vecteur des estimations
% et des erreurs. Pour trouver l'intervalle
% de confiance plus tard il suffirat d'uti-
% liser IC=[I-Z*s/sqrt(n),I+Z*s/sqrt(n)]
I_put=zeros(1,Nsim);
Err_put=I_put;

I_call=zeros(1,Nsim);
Err_call=I_call;

% parametres des options
K=1;
beta=1;

% les valeurs exactes des options
C_ex=(exp(beta^2/2)*normcdf(beta-log(K)/beta,0,1)-...
    K*normcdf(-log(K)/beta,0,1));
C=C_ex*ones(1,Nsim);

P_ex=(-exp(beta^2/2)*normcdf(-beta+log(K)/beta,0,1)+...
    K*normcdf(log(K)/beta,0,1));
P=P_ex*ones(1,Nsim);

% on augmente le nombre d'échantillons de pas à chaque
% fois pour obtenir Nsim simulation de chaque option

% debut de la boucle for pour les Nsim estimations
for n=1:Nsim    
    % estimation du call
    [I_hat_call,err_std_call]=monteCarloCall(n*pas);
    I_call(n)=I_hat_call;
    Err_call(n)=err_std_call;
    
    % estimation du put
    [I_hat_put,err_std_put]=monteCarloPut(n*pas);
    I_put(n)=I_hat_put;
    Err_put(n)=err_std_put;
end
%fin de la bourcle for pour les Nsim estimations

% affichage des resultats
% preparation de la figure
figure()
title('Monte carlo Exercice 1: question 1 2 et 3 avec \alpha =0.95')
xlabel('Nombre de simulations')
ylabel('Valeurs des options')
hold on

% affichage du put
figEx_put=plot(P);% valeur exacte
figI_put=plot(I_put);% l'estimation
figHaut_put=plot(I_put+Z*Err_put);% borne sup de l'IC
figBas_put=plot(I_put-Z*Err_put);% borne inf de l'IC

% affichage du call 
figEx_call=plot(C);% valeur exacte
figI_call=plot(I_call);% l'estimation
figHaut_call=plot(I_call+Z*Err_call);% borne sup de l'IC
figBas_call=plot(I_call-Z*Err_call);% borne inf de l'IC

% ajout des legendes à la figure
legend([figEx_put, figI_put, figHaut_put, figBas_put,...
        figEx_call, figI_call, figHaut_call, figBas_call],...
        'exact put', 'estimation put','haute put','basse call',...
        'exact call', 'estimation call','haute call','basse call');



%% test distance de Levenshtein

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

%% test de text sur une figure
A = 1/eps;
str_e = sprintf('ggas 18 %0.5e',A)