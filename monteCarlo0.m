function [I_hat,flag,err_std,nb_sim]=monteCarlo0(alpha,delta,N,g,f)
  % Cette fonction fait une simulation de Monte Carlo avec un niveau de
  % confiance alpha et une erreur absolue delta donnée également.
  % ENTREE : alpha: Niveau de confiance
  %          delta: Erreur absolue
  %          N: Le nombre maximum de simulation
  %          g: Fonction dont on veut calculer l'intégrale
  %          f: Densité de X
  % SORTIE : I_hat: La valeur approchée de la quantité que nous voulons
  %                 évaluer.
  %          flag: Variable "booleen" valant 1 si nous avons la precison 
  %                souhaitée sinon 0
  %          err_std: Erreur standard de la simulation réalisée
  %          nb_sim:  Nombre de simulations effectuées.
  
  X= pi*rand(); % nous supposons que nous evaluons une intégrande sur [0,pi]
  
  S1=X; % somme partielle des Yi
  
  S2=X^2; % somme patielle des carrées des Yi
  
  s=inf;% on intialise la variance avec l'infinie
  
  Z=norminv((alpha+1)/2,0,1); % le quantile d'ordre (alpha+1)/2 de la loi
  % normale
  
  % On fait un nombre n de simulation et on teste si la precision souhaitée
  % est atteinte.
  nb_sim=1;
 while(nb_sim<N && ((s/sqrt(nb_sim))*Z)>delta)
     
      nb_sim=nb_sim+1;
      
      X=pi*rand();% on simule une variable aléatoire uniforme sur [0,pi]
      Y=g(X)/f(X);% on évalue ensuite Y 
      S1=S1+Y;    % mise-à-jour de S1
      S2=S2+Y^2;  % mise-à-jour de S2
      %s=(S2-(S1/nb_sim)^2)/(nb_sim-1);% mise-à-jour de s
      s=(S2-nb_sim*(S1/nb_sim)^2)/(nb_sim-1);
      
 end
  
 flag=((s/sqrt(nb_sim))*Z)<=delta; % le flag c'est pour preciser si nous avons la
 % precision voulue
 
 % on retourne l'estimation obtenue.
 I_hat=S1/nb_sim;
 
 % on retourne l'erreur standard de cette simulation
 err_std=s/sqrt(nb_sim);
end
