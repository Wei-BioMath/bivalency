function run_ODE
% Output: steady-state solutions and phase plane for the deterministic model
%
% By Wei Zhao & Lingxia Qiao (2021)
k=[2 2];d=1*[1 1];N=60;K=[1 1]/30;gamma=0.2*[1 1];
nsample=200;
x0_temp=lhsu([0 0],[N N],nsample*2+100);x0=[];
for j=1:nsample*2+100                                                                                                                                                            
    if sum(x0_temp(j,:))<=N
        x0=[x0;[x0_temp(j,1) 0 0 0 0 x0_temp(j,2)]];
    end
end
x0=x0(1:nsample,:);x0(end,[1 6])=[0 60];x0(end-1,[1 6])=[60 0];x0(end-2,[1 6])=[30 30];x0(end-3,[1 6])=[0 0];
[state_stable,state_all,grad_absmax,lamda]=phase(x0,k,K,d,gamma,N);
end

function [state,state_all,grad_absmax,lamda]=phase(x0,k,K,d,gamma,N)
opt=optimset('Display','off');
%% random initials
AA=zeros(length(x0),6);
AA_all=AA;
figure;
for j=1:size(x0,1)
    T=1e5;
    [t,y]=ode15s(@(t,x)f(x,t,k,K,d,gamma,N),[0 T],x0(j,:));
    AA(j,1:6)=y(end,1:6);
    h1=plot(y(:,6),y(:,1),'-','color',[0.5 0.5 0.5]);
    hold on;
    h2=plot(y(end,6),y(end,1),'.r','MarkerSize',50);
    l=size(y,1);ll=round(l/5);
    % quiver to plot the direction of trajectory
    df=ff(y(ll,:),k,K,d,gamma,N);
    q=quiver(y(ll,6),y(ll,1),df(6),df(1),0.1);
    q.Color=[0.5 0.5 0.5];q.MaxHeadSize=1.5;q.LineWidth=1;
    [y,~,exitflag]=fsolve(@(x)ff(x,k,K,d,gamma,N),x0(j,:),opt);
    if exitflag==1
       plot(y(6),y(1),'or','MarkerSize',15);
       axis([-0.2 N+.2 -0.2 N+.2]);xlabel('H3K4me3');ylabel('H3K27me3');
       AA_all(j,1:6)=y;
    end
end
axis([0 60 0 60]);
xlabel('# of H3K4me3');ylabel('# of H3K27me3');
state=AA(1,:);
state_all=AA_all(1,:);
for i=2:length(AA)
    s=0;
    for j=1:size(state,1)
        if max(abs(AA(i,:)-state(j,:)))<0.1
            s=1;break;
        end
    end
    if s==0
       state=[state;AA(i,:)];
    end
end
for i=2:length(AA_all)
    s=0;
    for j=1:size(state_all,1)
        if max(abs(AA_all(i,:)-state_all(j,:)))<0.1
            s=1;break;
        end
    end
    if s==0 & min(AA_all(i,:))>0
       state_all=[state_all;AA_all(i,:)];
    end
end
grad_absmax=zeros(size(state_all,1),1);lamda=zeros(size(state_all,1),6);
% for j=1:size(state_all,1)
%     [lamdaj,gradj]=get_lamda(state_all(j,:),k,K,d,gamma,N);
%     grad_absmax(j)=gradj;
%     lamda(j,1:6)=lamdaj;
% end
end

function dy=f(x,t,k,K,d,gamma,N)
R3=x(1);R2=x(2);R1=x(3);
A1=x(4);A2=x(5);A3=x(6);
U=N-sum(x);
E=A3/N/(A3/N+K(1))+gamma(1);F=R3/N/(R3/N+K(2))+gamma(2);
v1=-d(2)*R3*E+k(2)*R2*F;
v2=-d(2)*R2*E+k(2)*R1*F;
v3=-d(2)*R1*E+k(2)*U*F;
v4=-k(1)*U*E+d(1)*A1*F;
v5=-k(1)*A1*E+d(1)*A2*F;
v6=-k(1)*A2*E+d(1)*A3*F;
dy(1)=v1;
dy(2)=-v1+v2;
dy(3)=-v2+v3;
dy(4)=-v4+v5;
dy(5)=-v5+v6;
dy(6)=-v6;
dy=dy';
end

function dy=ff(x,k,K,d,gamma,N)
R3=x(1);R2=x(2);R1=x(3);
A1=x(4);A2=x(5);A3=x(6);
U=N-sum(x);
E=A3/N/(A3/N+K(1))+gamma(1);F=R3/N/(R3/N+K(2))+gamma(2);
v1=-d(2)*R3*E+k(2)*R2*F;
v2=-d(2)*R2*E+k(2)*R1*F;
v3=-d(2)*R1*E+k(2)*U*F;
v4=-k(1)*U*E+d(1)*A1*F;
v5=-k(1)*A1*E+d(1)*A2*F;
v6=-k(1)*A2*E+d(1)*A3*F;
dy(1)=v1;
dy(2)=-v1+v2;
dy(3)=-v2+v3;
dy(4)=-v4+v5;
dy(5)=-v5+v6;
dy(6)=-v6;
dy=dy';
end

function [lamda,grad]=get_lamda(x,k,K,d,gamma,N)
syms R3 R2 R1 A1 A2 A3;
U=N-sum([R3 R2 R1 A1 A2 A3]);
E=A3/N/(A3/N+K(1))+gamma(1);F=R3/N/(R3/N+K(2))+gamma(2);
v1=-d(2)*R3*E+k(2)*R2*F;
v2=-d(2)*R2*E+k(2)*R1*F;
v3=-d(2)*R1*E+k(2)*U*F;
v4=-k(1)*U*E+d(1)*A1*F;
v5=-k(1)*A1*E+d(1)*A2*F;
v6=-k(1)*A2*E+d(1)*A3*F;
dy(1)=v1;
dy(2)=-v1+v2;
dy(3)=-v2+v3;
dy(4)=-v4+v5;
dy(5)=-v5+v6;
dy(6)=-v6;
dy=dy';
J=jacobian(dy,[R3 R2 R1 A1 A2 A3]);
R3=x(1);R2=x(2);R1=x(3);A1=x(4);A2=x(5);A3=x(6);
grad=double(max(abs(subs(dy))));
[V,D]=eig(subs(J));
lamda=double(diag(D))';
end

function s=lhsu(xmin,xmax,nsample)
% s=lhsu(xmin,xmax,nsample)
% LHS from uniform distribution
% Input:
%   xmin    : min of data (1,nvar)
%   xmax    : max of data (1,nvar)
%   nsample : no. of samples
% Output:
%   s       : random sample (nsample,nvar)
%   Budiman (2003)
rng('shuffle');
nvar=length(xmin);
ran=rand(nsample,nvar);
s=zeros(nsample,nvar);
for j=1: nvar
   idx=randperm(nsample);
   P =(idx'-ran(:,j))/nsample;
   s(:,j) = xmin(j) + P.* (xmax(j)-xmin(j));
end
end