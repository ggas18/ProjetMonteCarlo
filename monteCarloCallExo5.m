function [I_hat,err_std]=monteCarloCallExo5(N)
  % Cette fonction fait une simulation de Monte en
  % N echantillons en utilisant la symetrie de la
  % loi normale centree reduite pour faire des v.a
  % antithetiques
  % ENTREE : N: Le nombre de simulation      
  % SORTIE : I_hat: La valeur approchee du call 
  %          err_std: Erreur standard de la simulation realisee
  %                   parametre du call
  beta=1;K=1;
  X= randn(); % simulation d'une variable de loi normale
  % centree reduite
  Y_1=max(exp(beta*X)-K,0);% antithetique 1
  Y_2=max(exp(-beta*X)-K,0);% antithetique 2
  S1_1=Y_1; % somme partielle des Y1i
  S2_1=Y_1^2; % somme patielle des Y1i^2
  S1_2=Y_2; % somme partielle des Y2i
  S2_2=Y_2^2; % somme patielle des Y2i^2
  m_y1y2=Y_1*Y_2;% la somme des produits des Y1i*Y2i
  % utile pour calculer la covariance
  n=1;
 while(n<N)
      X=randn();% on une normale centree reduite
      Y_1=max(exp(beta*X)-K,0);% antithetique 1
      Y_2=max(exp(-beta*X)-K,0);% antithetique 2
      S1_1=S1_1+Y_1;    % mise-a-jour de S1_1
      S2_1=S2_1+Y_1^2;  % mise-a-jour de S2_1 
      S1_2=S1_2+Y_2;    % mise-a-jour de S1_2
      S2_2=S2_2+Y_2^2;  % mise-a-jour de S2_2 
      m_y1y2= m_y1y2+Y_1*Y_2;% mise-a-jour de m_y1y2
      n=n+1; 
 end
 % on estime la variance par son estimateur sans biais
 % pour la v.a Y1
 v_1=(S2_1-N*(S1_1/N)^2)/(N-1);
 % on estime la variance par son estimateur sans biais
 % pour la v.a Y2
 v_2=(S2_2-N*(S1_2/N)^2)/(N-1);
 % la covariance
 cov1_2=(m_y1y2-(S1_1*S1_2)/N)/(N-1);
 % variance et ecart-type finaux est
 variance=v_1/4+v_2/4+cov1_2/2;
 s=sqrt(variance);
 
 % on retourne l'estimation obtenue.
 I_hat=(S1_1+S1_2)/(2*N);
 % on retourne l'erreur standard de cette simulation
 err_std=s/sqrt(N);
end