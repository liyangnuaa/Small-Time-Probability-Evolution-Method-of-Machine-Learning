clear;
clc;

xmin=-1;
xmax=1;
zmin=-5;
zmax=5;
Nz=5e6;
z0=linspace(zmin,zmax,Nz);

% Cinterval=1;             %%% 计算漂移系数的积分区间
xlim=100;                %%% 跳跃振幅截断
% N=1;                   %%% 一个初始点给出的样本数
h=0.001;
alpha=1.5;
sigma=2;

M=h^(1/alpha)*stblrnd(alpha,0,1,0,1,Nz);
xf=(4*z0-z0.^3)*h+sigma*M;
xf0=z0;

a=0.2;
m=5;
N=2;
nk=zeros(1,N+1);
for k=0:N
    I=(abs(xf)>=m^k*a)&(abs(xf)<m^(k+1)*a);
    nk(k+1)=length(xf(I));
end
nratio=nk(1)./nk(2:end);
pos=1:N;
alpha1=log(nratio)./(pos*log(m));
alpha0=sum(alpha1)/N;

cnalpha=alpha0*gamma((1+alpha0)/2)/(2^(1-alpha0)*sqrt(pi)*gamma(1-alpha0/2));
pos=0:N;
sigmak=a*m.^pos.*(nk*alpha0/(h*Nz*2*cnalpha*(1-m^(-alpha0)))).^(1/alpha0);
sigma0=sum(sigmak)/(N+1);

% xc=0.1:0.01:1.9;
% cnalpha=xc.*gamma((1+xc)/2)./(2.^(1-xc)*sqrt(pi).*gamma(1-xc/2));
% figure;
% plot(xc,cnalpha);
% 
Ncoef=6;
I=(abs(xf)<xlim);
zinitial=xf0(I);
x=xf(I);
n=length(zinitial);
A=zeros(n,Ncoef+1);
A(:,1)=1;
for i=1:Ncoef
    A(:,i+1)=zinitial'.^i;
end
B=x'/h;
% X=(A'*A)\(A'*B);
pos=0:1:Ncoef;
for k=1:Ncoef
    X=(A'*A)\(A'*B);
    I=abs(X)<0.1;
    A(:,I)=[];
    pos(I)=[];
    if isempty(X(I))
        break;
    end
end

