%% monte carlo pour l'exercice 2
clc;close all
% le quantile d'ordre (alpha+1)/2 de la loi normale
% centree reduite.
alpha=.95;
Z=norminv((alpha+1)/2,0,1); 
% pas et nombre d'estimations des options.
Nsim=200;pas=200;
% on initialise le vecteur des estimations
% et des erreurs. Pour trouver l'intervalle
% de confiance plus tard il suffirat d'uti-
% liser IC=[I-Z*s/sqrt(n),I+Z*s/sqrt(n)]
I_put=zeros(1,Nsim);
Err_put=I_put;
I_call_1=zeros(1,Nsim);
Err_call_1=I_call_1;
I_call_2=zeros(1,Nsim);
Err_call_2=I_call_2;
% parametres des options
K=1;beta=1;
% les valeurs exactes des options
C_ex=(exp(beta^2/2)*normcdf(beta-log(K)/beta,0,1)-...
    K*normcdf(-log(K)/beta,0,1));
C=C_ex*ones(1,Nsim);
% on augmente le nombre d'echantillons de pas a chaque
% fois pour obtenir Nsim simulation de chaque option

% debut de la boucle for pour les Nsim estimations
% estimation du call avec la methode de l'exo1
start=tic;
for n=1:Nsim    
    [I_hat_call_1,err_std_call_1]=monteCarloCall(n*pas);
    I_call_1(n)=I_hat_call_1;
    Err_call_1(n)=err_std_call_1;   
end
t_exo1=toc(start);
% fin de la boucle for pour les Nsim estimations exo1
% debut de la boucle for pour les Nsim estimations
% estimation du call avec la methode de l'exo2
start=tic;
for n=1:Nsim     
    [I_hat_call_2,err_std_call_2]=monteCarloCallExo2(n*pas);
    I_call_2(n)=I_hat_call_2;Err_call_2(n)=err_std_call_2;   
end
t_exo2=toc(start);
% fin de la boucle for pour les Nsim estimations exo1

% affichage des resultats et preparation de la figure
N=pas*(1:Nsim); % vecteur du nombre de simulation
fig_exo2=figure();
title('Exo 2 call seul: question 2 et 3 avec \alpha =0.95')
xlabel('Nombre de simulations')
ylabel('Valeurs des options')
hold on
% ajout des temps de calculs
str = {sprintf('temps exo 1: %0.2es',t_exo1),...
       sprintf('temps exo 2: %0.2es',t_exo2)};
text(pas*Nsim*7/14,.5,str)
% affichage du call exact
figEx_call=plot(N,C,'LineWidth',2);% valeur exacte
% affichages pour la methode de l'exo2
figI_call_2=plot(N,I_call_2,'LineWidth',2);% l'estimation
figHaut_call_2=plot(N,I_call_2+Z*Err_call_2,'LineWidth',2);% sup de l'IC
figBas_call_2=plot(N,I_call_2-Z*Err_call_2,'LineWidth',2);%inf de l'IC
% affichage du call par l'exo 1 
figI_call_1=plot(N,I_call_1,'LineWidth',2);% l'estimation
figHaut_call_1=plot(N,I_call_1+Z*Err_call_1,'LineWidth',2);% borne sup de l'IC
figBas_call_1=plot(N,I_call_1-Z*Err_call_1,'LineWidth',2);% borne inf de l'IC
% ajout des legendes a la figure
legend([ figI_call_2, figHaut_call_2, figBas_call_2,...
        figEx_call, figI_call_1, figHaut_call_1, figBas_call_1],...
         'estimation exo2','haute exo2','basse exo2',...
        'exact call', 'estimation exo1','haute exo1','basse exo1');
 %%
 % on enregistre la figure sous format jpg
chem='images';
chem=strcat(chem,'/exo2');
print(fig_exo2,chem,'-djpeg')
