% perform stochastic simulation algorithm (SSA) to obtain the stochastic trajectory  
clear all;
%% parameters
% see parameter meaning by type "help main" in command window or look at
% the main function directly
MAX_ITER=2e5;period=20;N=60;NS=7;
local=0;division=1;NSA=(NS-1)/2; % the number of (de)methylation
ka1=2;kr1=2;gamma=0.2;% gamma is the noise-to-feedback ratio
da=1*ones(NSA,1);eda=gamma*da;   
dr=1*ones(NSA,1);edr=gamma*dr;  
ka=ka1*ones(NSA,1);eka=gamma*ka; 
kr=kr1*ones(NSA,1);ekr=gamma*kr;
rho=0; 
%% initialization
X0=3*ones(N,1);
P=zeros(1,4); % time-averaged probability of the four states
nP=zeros(1,3); % time-averaged probability of nucleosomes with K4me3/K4me3, K4me3/K27me3 and K27me3/K27me3
%% perform the Gillespie algoritm
[P,nP]=main_SSA(@fnl,ka,kr,da,dr,eka,ekr,eda,edr,rho,...
    X0,MAX_ITER,period,N,NS,local,division);