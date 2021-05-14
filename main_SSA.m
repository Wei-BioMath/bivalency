function [P,nP]=main_SSA(f,ka,kr,da,dr,eka,ekr,eda,edr,rho,X,MAX_ITER,period,N,NS,local,division)
% Perform stochastic simulation algorithm (SSA) to obtain the stochastic trajectory
% 
% Output: time-averaged probability of four chromatin states and three
% kinds of nucleosomes
%
% Model settings:
% f: function handle @fnl means nonlinear feedback and @fl means linear feedback
% N: histone number 
% NS: number of H3 histone states. For 7-state model, set NS=7.
% local: loca=1 means local interaction regime, then nei=5; local=0 means global interaction regime, then nei=N-1
% division: division=1 means DNA replication is applied.
% period: cell cycle length.
% MAX_ITER: maximum iteractions for Gillespie algorithm
%
% Kinetic parameters:
% ka: feedback-mediated H3K4 methylation rate constants
% eka: noisy H3K4 methylation rates
% da: feedback-mediated H3K4 demethylation rate constants
% eda: noisy H3K4 demethylation rates
% kr: feedback-mediated H3K27 methylation rate constants
% ekr: noisy H3K27 methylation rates
% da: feedback-mediated H3K27 demethylation rate constants
% eda: noisy H3K27 demethylation ratess
%
% X: initial condition
% rho: reduced A2 or R2 mediated positive feedback strength, not used.
%
% By Wei Zhao & Lingxia Qiao (2021)

K=[2 2]/60; % the Hill coefficient used for the nonlinear feedback

tic
%% initialization
NSA=(NS-1)/2;
T=zeros(MAX_ITER,1); % store time points
XX=zeros(N,MAX_ITER); % store state vectors
EA3=zeros(N,1);EA2=zeros(N,1); % summing the positive feedback from active histone A3 or A2
EA=zeros(N,1); % store numbers of H3K4me3 histone neighbors for each histone
FR3=zeros(N,1);FR2=zeros(N,1); % summing the positive feedback from repressive histone R3 or R2
FR=zeros(N,1); % store numbers of H3K27me3 histone neighbors for each histone
fa=zeros(N,1);fr=zeros(N,1); % store nonlinear or linear feedback for each histone
A=zeros(MAX_ITER,NSA);% the number of A3 A2 A1, i.e., K4me3,K4me2,K4me1
R=zeros(MAX_ITER,NSA);% the number of R3 R2 R1, i.e., K27me3,K27me2,K27me1
iter=1;XX(:,1)=X;
A(iter,1)=length(find(X==NSA));A(iter,2)=length(find(X==NSA-1));A(iter,3)=length(find(X==NSA-2));
R(iter,1)=length(find(X==-NSA));R(iter,2)=length(find(X==-NSA+1));R(iter,3)=length(find(X==-NSA+2));
% E, F initialization
for j=2:2:N
    if local==1
        nei=5;neighbor=[j-3 j-2 j-1 j+1 j+2];neighbor(neighbor<=0 | neighbor>N)=[];
        elseif local==0
        nei=N-1;neighbor=1:N;neighbor(j)=[];
    end
    EA3(j)=sum(eq(X(neighbor),NSA));EA2(j)=sum(eq(X(neighbor),NSA-1));
    FR3(j)=sum(eq(X(neighbor),-NSA));FR2(j)=sum(eq(X(neighbor),-NSA+1));
end
for j=1:2:N
    if local==1
        neighbor=[j-2 j-1 j+1 j+2 j+3];neighbor(neighbor<=0 | neighbor>N)=[];
    elseif local==0
        nei=N-1;neighbor=1:N;neighbor(j)=[];
    end 
    EA3(j)=sum(eq(X(neighbor),NSA));EA2(j)=sum(eq(X(neighbor),NSA-1));
    FR3(j)=sum(eq(X(neighbor),-NSA));FR2(j)=sum(eq(X(neighbor),-NSA+1));
end
EA=EA3+EA2*rho;FR=FR3+FR2*rho;
fa=f(EA,FR,nei,K(1:length(K)/2));fr=f(FR,EA,nei,K(length(K)/2+1:end));
% rate/propensity initialization
p=fa*[dr;ka]'+ones(N,1)*[edr;eka]'; % pos_direction 
m=fr*[da;kr]'+ones(N,1)*[eda;ekr]'; % neg_direction
rate=zeros(N,2);
% eq(X,-NSA); p(:,1),
% eq(X,-NSA+1); p(:,2),m(:,end) 
% ... 
% eq(X,NSA-1);p(:,2*NSA), m(:,2) 
% eq(X,NSA); m(:,1)
for j=1:2*NSA
    rate(:,1)=rate(:,1)+eq(X,-NSA+j-1).*p(:,j);
    rate(:,2)=rate(:,2)+eq(X,-NSA+j).*m(:,end-j+1);
end
%% iteration
RN1=rand(MAX_ITER,1);
RN2=rand(MAX_ITER,1); % random numbers 
while iter<MAX_ITER
r_long=rate';r_long=r_long(:);r_tot=sum(r_long);rn1=RN1(iter);rn2=RN2(iter);
t=-log(rn1)/r_tot;
if division && floor(T(iter)/period)<floor((T(iter)+t)/period) % if division then true; false otherwise
    % binomial partition
    RN3=rand(N/2,1)<0.5;RN4=reshape([RN3 RN3]',[N,1]);
    X(RN4==0)=0;iter=iter+1;
    XX(:,iter)=X;T(iter)=floor((T(iter-1)+t)/period)*period;
    A(iter,1)=length(find(X==NSA));A(iter,2)=length(find(X==NSA-1));A(iter,3)=length(find(X==NSA-2));
    R(iter,1)=length(find(X==-NSA));R(iter,2)=length(find(X==-NSA+1));R(iter,3)=length(find(X==-NSA+2));
    for j=2:2:N
        if local==1
            neighbor=[j-3 j-2 j-1 j+1 j+2];neighbor(neighbor<=0 | neighbor>N)=[];
            elseif local==0
            neighbor=1:N;neighbor(j)=[];
        end
        EA3(j)=sum(eq(X(neighbor),NSA));EA2(j)=sum(eq(X(neighbor),NSA-1));
        FR3(j)=sum(eq(X(neighbor),-NSA));FR2(j)=sum(eq(X(neighbor),-NSA+1));
    end
    for j=1:2:N
        if local==1
            neighbor=[j-2 j-1 j+1 j+2 j+3];neighbor(neighbor<=0 | neighbor>N)=[];
        elseif local==0
            nei=N-1;neighbor=1:N;neighbor(j)=[];
        end 
        EA3(j)=sum(eq(X(neighbor),NSA));EA2(j)=sum(eq(X(neighbor),NSA-1));
        FR3(j)=sum(eq(X(neighbor),-NSA));FR2(j)=sum(eq(X(neighbor),-NSA+1));
    end
    EA=EA3+EA2*rho;FR=FR3+FR2*rho;fa=f(EA,FR,nei,K(1:length(K)/2));fr=f(FR,EA,nei,K(length(K)/2+1:end));
    % rate/propensity initialization
    p=fa*[dr;ka]'+ones(N,1)*[edr;eka]'; % pos_direction 
    m=fr*[da;kr]'+ones(N,1)*[eda;ekr]'; % neg_direction
    rate=zeros(N,2);
    for j=1:2*NSA
        rate(:,1)=rate(:,1)+eq(X,-NSA+j-1).*p(:,j);
        rate(:,2)=rate(:,2)+eq(X,-NSA+j).*m(:,end-j+1);
    end
else
    r_cum=cumsum(r_long);
    loc=find((rn2*r_tot-r_cum)<0);loc=loc(1);
    ith=floor((loc-1)/2)+1; % update ith: 1:N
    jth=loc-(ith-1)*2;% update jth: 1 or 2 ; 1 means +1 reaction,2 means -1 reaction
    X(ith)=X(ith)-2*jth+3;iter=iter+1;
    XX(:,iter)=X;T(iter)=T(iter-1)+t;
    A(iter,1)=length(find(X==NSA));A(iter,2)=length(find(X==NSA-1));A(iter,3)=length(find(X==NSA-2));
    R(iter,1)=length(find(X==-NSA));R(iter,2)=length(find(X==-NSA+1));R(iter,3)=length(find(X==-NSA+2));
    % E, F update & rate propensity update
    if local==1
        if mod(ith,2)==1
        neighbor=[ith-2 ith-1 ith+1 ith+2 ith+3];
        else 
        neighbor=[ith-3 ith-2 ith-1 ith+1 ith+2];
        end
        neighbor(neighbor<=0 | neighbor>N)=[];
    elseif local==0
        neighbor=1:N;neighbor(ith)=[];
    end
    EA(neighbor)=EA(neighbor)...
        +eq(X(ith),NSA)*(1-rho)...
        +eq(X(ith),NSA-1)*eq(jth,2)*(rho-1)...
        +eq(X(ith),NSA-1)*eq(jth,1)*rho...
        -eq(X(ith),NSA-2)*eq(jth,2)*rho;
    FR(neighbor)=FR(neighbor)...
    +eq(X(ith),-NSA)*(1-rho)...
    +eq(X(ith),-NSA+1)*eq(jth,1)*(rho-1)...
    +eq(X(ith),-NSA+1)*eq(jth,2)*rho...
    -eq(X(ith),-NSA+2)*eq(jth,1)*rho;
    % rate/propensity update
    p(neighbor,:)=f(EA(neighbor),FR(neighbor),nei,K(1:length(K)/2))*[dr;ka]'+ones(length(neighbor),1)*[edr;eka]'; % pos_direction 
    m(neighbor,:)=f(FR(neighbor),EA(neighbor),nei,K(length(K)/2+1:end))*[da;kr]'+ones(length(neighbor),1)*[eda;ekr]'; % neg_direction
    rate([neighbor ith],1:2)=0;
    for j=1:2*NSA
        rate([neighbor ith],1)=rate([neighbor ith],1)+eq(X([neighbor ith]),-NSA+j-1).*p([neighbor ith],j);
        rate([neighbor ith],2)=rate([neighbor ith],2)+eq(X([neighbor ith]),-NSA+j).*m([neighbor ith],end-j+1);
    end
end
end
%% time-averaged probability Pa, Pr, Pb, Pu
A_trim=A(1:end-1,1);
R_trim=R(1:end-1,1);
T_trim=T(1:end-1);
T_weigth=T(2:end)-T(1:end-1);
PA=T_weigth'*(A_trim(:,1)>2*R_trim(:,1) & A_trim(:,1).^2+R_trim(:,1).^2>(N/4)^2)/(T_trim(end)-T_trim(1));
PR=T_weigth'*(R_trim(:,1)>2*A_trim(:,1) & A_trim(:,1).^2+R_trim(:,1).^2>(N/4)^2)/(T_trim(end)-T_trim(1));
PB=T_weigth'*(A_trim(:,1).^2+R_trim(:,1).^2>(N/4).^2 & A_trim(:,1)<=2*R_trim(:,1) & R_trim(:,1)<=2*A_trim(:,1))/(T_trim(end)-T_trim(1));
PU=T_weigth'*(A_trim(:,1).^2+R_trim(:,1).^2<=(N/4).^2)/(T_trim(end)-T_trim(1));
P=[PA PR PB PU];

%% time-averaged probability of nucleosomes at the bivalent state
XX_trim=XX(:,1:end-1);
ind_biv=find(A_trim(:,1).^2+R_trim(:,1).^2>(N/4)^2 ...
    & A_trim(:,1)<=2*R_trim(:,1) & R_trim(:,1)<=2*A_trim(:,1));% searching the bivalent state first
nP=zeros(1,3);% probability of nucleosomes with K4me3/K4me3, K4me3/K27me3 and K27me3/K27me3
if ~isempty(ind_biv) 
XX_biv=XX_trim(:,ind_biv);
nPP=zeros(length(ind_biv),3); % number of of nucleosomes with K4me3/K4me3, K4me3/K27me3 and K27me3/K27me3 at each time point
for jjj=1:length(ind_biv)
    xt=XX_biv(:,jjj);
    xt=reshape(xt,2,N/2)';
    nPP(jjj,1)=sum((xt(:,1)==NSA).*(xt(:,2)==NSA));
    nPP(jjj,2)=sum((xt(:,1)==NSA).*(xt(:,2)==-NSA))+sum((xt(:,1)==-NSA).*(xt(:,2)==NSA));
    nPP(jjj,3)=sum((xt(:,1)==-NSA).*(xt(:,2)==-NSA));
end
nP=T_weigth(ind_biv)'*nPP/sum(T_weigth(ind_biv))/N*2*100;
end


%% plot the dynamics
figure;
XX_interp=interp1(T,XX',0:0.1:T(end),'previous');
subplot(121);im=imagesc(T/period,1:N,XX_interp');set(gca,'YDir','normal');
c=colorbar('Ticks',[-3:3],'TickLabels',{'H3K27me3','H3K27me2','H3K27me1','Empty','H3K4me1','H3K4me2','H3K4me3'});
ylabel("Histone position"); xlabel("Time (cell cycles)");
subplot(122);hold on;h1=plot(T/period,A(:,1),'-r');h2=plot(T/period,R(:,1),'-b');
xlabel("Time (cell cycles)");ylabel("Histone number");ylim([-0.2 N+0.2]);
legend([h1 h2],{'H3K4me2/3','H3K27me2/3'});

toc
end

function y=fnl(x,x1,N,K)
x=x/N;
y=x./(x+K);
end

function y=fl(x,x1,N,K)
y=x/N;
end
