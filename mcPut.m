function [I_hat,IC_h,IC_b,Op_exact]=mcPut(pas,N_t,alpha)
    % cette fonction prend deux entiers l'un est le pas et le
    % le second est le nombre d'estimations

    % ENTREE: pas: le pas des nombres de simulations
    %         N_t: le nombre d'estimations de la valeur
    %              que nous voulons approchées
    %         alpha: le niveau de l'intervalle de confiance
    %                ce n'est pas obligatoire de passer ce
    %                parametre, par defaut il est 0.95

    %
    % SORTIE: I_hat: le vecteur de taille N_t donnant les
    %                les valeurs que nous avons estimées
    %         IC_h:  le vecteur de taille N_t donnant les
    %                les bornes hautes de l'invertalle de 
    %                confiance
    %         IC_b:  le vecteur de taille N_t donnant les
    %                les bornes basses de l'invertalle de 
    %                confiance
    %         C_exact:  un vecteur dont chaque composante
    %                vaut la valeur exacte de l'options
    %                
    
    % les paramètres de l'option
    beta=1;
    K=1;
    
    % si alpha n'est pas donné alors on fait alpha=0.95
    if(nargin <3) 
        alpha=0.95;
    end
    
    Op_exact=(-exp(beta^2/2)*normcdf(-beta+log(K)/beta,0,1)+K*normcdf(log(K)/beta,0,1))*ones(1,N_t);
    
    % le quantile d'ordre (alpha+1)/2 de la loi normale
    % centrée reduite.
    Z=norminv((alpha+1)/2,0,1); 
    
    %on initialise les vecteurs de retour à 0
    I_hat=zeros(1,N_t);
    IC_h=zeros(1,N_t);
    IC_b=zeros(1,N_t);
    
    % simulation de pas variables normales centrées
    % reduites
    X=randn(1,pas);
    
    % on calcule le vecteur des valeurs de l'option en
    % vectorisant les calculs
    OP_sim=max(K-exp(beta*X),0);
    
    % a chaque pas, on calcule la moyenne et la variance
    S1=sum(OP_sim);
    S2=sum(OP_sim.^2);
    n=pas;
    for i=1:N_t
        % simulation de pas variables normales centrées
        % reduites
        X=randn(1,pas);
    
        % on calcule le vecteur des valeurs de l'option en
        % vectorisant les calculs
        OP_sim=max(K-exp(beta*X),0);
        
        
        I_hat(i)=S1/n;
        s=(S2-n*(I_hat(i))^2)/(n-1);
        IC_b(i)=I_hat(i)-Z*s/sqrt(n);
        IC_h(i)=I_hat(i)+Z*s/sqrt(n);
        
        n=n+pas;
        S1=S1+sum(OP_sim);
        S2=S2+sum(OP_sim.^2);
        
        
    end
end


