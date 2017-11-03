%% monte carlo pour l'exercice 1

clc
close all

% le quantile d'ordre (alpha+1)/2 de la loi normale
% centree reduite.
alpha=.95;
Z=norminv((alpha+1)/2,0,1); 

% le nombre d'estimations des options.
Nsim=200;
% les pas entre la taille des echantillons
pas=100;

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

% on augmente le nombre d'echantillons de pas a chaque
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
N=pas*(1:Nsim); % vecteur du nombre de simulation
fig_exo1=figure();
title('Monte carlo Exercice 1: question 1 2 et 3 avec \alpha =0.95')
xlabel('Nombre de simulations')
ylabel('Valeurs des options')
hold on

% affichage du put
figEx_put=plot(N,P,'LineWidth',.74);% valeur exacte
figI_put=plot(N,I_put,'LineWidth',.74);% l'estimation
figHaut_put=plot(N,I_put+Z*Err_put,'LineWidth',.74);% borne sup de l'IC
figBas_put=plot(N,I_put-Z*Err_put,'LineWidth',.74);% borne inf de l'IC

% affichage du call 
figEx_call=plot(N,C,'LineWidth',2);% valeur exacte
figI_call=plot(N,I_call,'LineWidth',2);% l'estimation
figHaut_call=plot(N,I_call+Z*Err_call,'LineWidth',2);% borne sup de l'IC
figBas_call=plot(N,I_call-Z*Err_call,'LineWidth',2);% borne inf de l'IC

% ajout des legendes a la figure
legend([figEx_put, figI_put, figHaut_put, figBas_put,...
        figEx_call, figI_call, figHaut_call, figBas_call],...
        'exact put', 'estimation put','haute put','basse call',...
        'exact call', 'estimation call','haute call','basse call');
 
 %%
 % on enregistre la figure sous format jpg
chem='images';
chem=strcat(chem,'/exo1');
print(fig_exo1,chem,'-djpeg')
