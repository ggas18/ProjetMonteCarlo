function [I_hat,err_std]=monteCarloCallExo2(N)
  % Cette fonction fait une simulation de Monte en
  % N echantillons en utilisant la deuxième expression
  % du call et on simule une loi exponentielle
  % ENTREE : N: Le nombre de simulation
  %          
  % SORTIE : I_hat: La valeur approchee de la quantite 
  %                  que nous voulons evaluer.
  %          err_std: Erreur standard de la simulation realisee
  %                   parametre du call
  
  X= exprnd(1); % simulation d'une variable loi exponentielle
  % de parametre 1
  Y=(exp(sqrt(2*X))-1)/sqrt(2*X);
  S1=Y; % somme partielle des Yi
  
  S2=Y^2; % somme patielle des Yi^2
  
  n=1;
  
 while(n<N)
      
      X=exprnd(1);% on une variable de loi exponentielle de parametre 1
      Y=(exp(sqrt(2*X))-1)/(sqrt(2*X)*sqrt(2*pi));% on evalue ensuite Y 
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
