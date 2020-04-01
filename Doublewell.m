clear;
clc;

xmin=-1;
xmax=1;
zmin=-5;
zmax=5;
Nz=1e8;
z0=linspace(zmin,zmax,Nz);

% Cinterval=1;             %%% 计算漂移系数的积分区间
xlim=100;                %%% 跳跃振幅截断
N=1;                   %%% 一个初始点给出的样本数
h=0.001;
alpha=0.5;
sigma=2;
xf0=zeros(1,Nz*N);
xf=zeros(1,Nz*N);
position=(0:Nz-1)*N;

M=h^(1/alpha)*stblrnd(alpha,0,1,0,1,Nz+1e4);
I=abs(M)<xlim;
M=M(I);
M=M(1:Nz);
xf(position+1)=(4*z0-z0.^3)*h+sigma*M;
xf0(position+1)=z0;

xpts=linspace(-xlim,xlim,1e3);
[pdf,x00]=ksdensity(xf,xpts);
pdf1=pdf/h;
x2=x00;
pdf2=pdf1;
[~,I]=min(abs(x00));
x2(I)=[];
pdf2(I)=[];
[~,I]=min(abs(x00));
x2(I)=[];
pdf2(I)=[];

syms t
% f=fittype('k1/(t+r1)^r2+k2','independent','t','coefficients',{'r1','r2','k1','k2'});
% cfun=fit(x0,y0,f,'StartPoint',[3,4,-1e5,38]); %显示拟合函数，数据必须为列向量形式
f=fittype('calpha/abs(t).^(1+alpha0)','independent','t','coefficients',{'calpha','alpha0'});
cfun=fit(x2',pdf2',f,'StartPoint',[0.3,1]); %显示拟合函数，数据必须为列向量形式
delta=0.05;
xi1=-xlim:0.01:-delta;
xi2=delta:0.01:xlim;
xi=[xi1 xi2];
% xi=-xlim:0.01:xlim+0.001;
yi=cfun(xi);
figure;
plot(x00,pdf1,'r*',xi,yi,'b-');

xc=0.1:0.01:1.9;
cnalpha=xc.*gamma((1+xc)/2)./(2.^(1-xc)*sqrt(pi).*gamma(1-xc/2));
figure;
plot(xc,cnalpha);

Ncoef=6;
I=(xf0<zmax)&(xf0>zmin);
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

