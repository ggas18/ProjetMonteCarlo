function [I_hat,err_std]=monteCarloCall(N)
  % Cette fonction fait une simulation de Monte en
  % N echantillons
  % ENTREE : N: Le nombre de simulation
  %          
  % SORTIE : I_hat: La valeur approchee de la quantite 
  %                  que nous voulons evaluer.
  %          err_std: Erreur standard de la simulation realisee
  %                   parametre du call
  beta=1;
  K=1;
  
  X= randn(); % simulation d'une variable de loi normale
  % centree reduite
  Y=max(exp(beta*X)-K,0);% on evalue ensuite Y 
  S1=Y; % somme partielle des Yi
  
  S2=Y^2; % somme patielle des Yi^2
  
  n=1;
  
 while(n<N)
      
      X=randn();% on simule une variable aleatoire uniforme sur [0,pi]
      Y=max(exp(beta*X)-K,0);% on evalue ensuite Y 
      S1=S1+Y;    % mise-a-jour de S1
      S2=S2+Y^2;  % mise-a-jour de S2
      
      n=n+1;
      
 end
 % on estime la variance par son estimateur sans biais
 s=sqrt((S2-N*(S1/N)^2)/(N-1));
 % on retourne l'estimation obtenue.
 I_hat=S1/N;
 % on retourne l'erreur standard de cette simulation
 err_std=s/sqrt(N);
end
