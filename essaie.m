%% Pour tester la génération de loi uniforme sur {0,..,n-1} en utilisant la méthode générale
n=100;
n_sim=2000;

X=zeros(n_sim,1);

for i=1:n_sim
    X(i)=uniformGen(n);
end

moy=mean(X);
var=std(X)^2;

%% Test des taux d'intérêt
clc
r=0:0.001:20/100;
t=[(1+r).^4;(1+r).^3; (1+r).^2; (1+r).^1; (1+r).^0];
Ax=[13; 14; 16; 18; 20];
B=[16; 16; 16; 16; 16];
C=[17; 17; 16; 16; 13];
Comp=[Ax'*t; B'*t; C'*t];

%%
i=3.7/100;
a1=10000;
r=a1/10;
val=((a1+r/i) *((i+1)^n - 1 )/i - n*r/i)*(1+i);

%% test de la fonction monteCarlo
clc
clear

g=@(x)log(x);
f=@(x) 1/pi;
[I_hat,flag,err_std,nb_sim]=monteCarlo(0.95,0.01,10000000,g,f);

%% Theorie des jeux: jeux à somme nulle
Aj=[4 2
    1 3];
f=[ones(1,size(Aj,1)) 1];

Ax=[-Aj  ones(size(Aj,1),1)];
Ay=[-Aj' ones(size(Aj,1),1)];
Aeq=[ones(1,size(Aj,1)) 0];
beq=1;

b=zeros(1,size(Aj,2));
x=linprog(-f,Ax,b,Aeq,beq);
y=linprog(-f,Ay,b,Aeq,beq);
