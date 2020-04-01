clear;
clc;

zm=2;
zxmin=-zm;
zxmax=zm;
zymin=-zm;
zymax=zm;
Nzx=1e4;
Nzy=1e4;
Nz=Nzx*Nzy;
y0=linspace(zymin,zymax,Nzy);
z0x=zeros(1,Nz);
z0y=zeros(1,Nz);
for i=1:Nzy
    z0x((i-1)*Nzx+1:i*Nzx)=linspace(zxmin,zxmax,Nzx);
    z0y((i-1)*Nzx+1:i*Nzx)=y0(i);
end

% Cinterval=1;             %%% 计算漂移系数的积分区间
xlim=100;                %%% 跳跃振幅截断
N=1;                   %%% 一个初始点给出的样本数
h=0.001;
alpha=0.5;
sigma=0.5;

% zxf0=z0x;
% zyf0=z0y;
zxf=zeros(1,Nz*N);
zyf=zeros(1,Nz*N);
zxf0=zeros(1,Nz*N);
zyf0=zeros(1,Nz*N);
position=(0:Nz-1)*N;

gamma0=gamma((1+alpha)/2)*gamma(1)*h/(sqrt(pi)*gamma(1+alpha/2));
M=stblrnd(alpha/2,1,2*(gamma0*cos(pi*alpha/4))^(2/alpha),0,1,Nz+1e4);
I=abs(M)<xlim^2;
M=M(I);
M=M(1:Nz);
Normal=randn(2,Nz);
zxf(position+1)=(z0x-z0x.^3-5*z0x.*z0y.^2)*h+sigma*sqrt(M).*Normal(1,:);
zxf0(position+1)=z0x;
zyf(position+1)=-(1+z0x.^2).*z0y*h+sigma*sqrt(M).*Normal(2,:);
zyf0(position+1)=z0y;

% zxf0=zxf;
% zyf0=zyf;
znorm=sqrt(zxf.^2+zyf.^2);
xf=znorm;

% xpts=linspace(1e-1,1e1,1e3).^2;
% xpts=linspace(-xlim,xlim,1e4);
xpts=linspace(1e-2,xlim,1e3);
[pdf,x00]=ksdensity(xf,xpts);
% pdf1=2*pdf/(h*2*pi);
pdf1=pdf/(h*2*pi);
x2=x00;
pdf2=pdf1;
[~,I]=min(abs(x00));
% I=1:2;
x2(I)=[];
pdf2(I)=[];
% [~,I]=min(abs(x00));
% x2(I)=[];
% pdf2(I)=[];


syms t
% f=fittype('k1/(t+r1)^r2+k2','independent','t','coefficients',{'r1','r2','k1','k2'});
% cfun=fit(x0,y0,f,'StartPoint',[3,4,-1e5,38]); %显示拟合函数，数据必须为列向量形式
f=fittype('calpha/abs(t).^(1+alpha0)','independent','t','coefficients',{'calpha','alpha0'});
cfun=fit(x2',pdf2',f,'StartPoint',[0.15,1]); %显示拟合函数，数据必须为列向量形式
delta=0.05;
% xi1=-xlim:0.01:-delta;
xi2=delta:0.01:xlim;
% xi=[xi1 xi2];
xi=xi2;
% xi=-xlim:0.01:xlim+0.001;
yi=cfun(xi);
figure;
plot(x00,pdf1,'r*',xi,yi,'b-');

% figure;
% plot(zxf0(1:1e7),zyf0(1:1e7),'.');
% axis([-xlim xlim -xlim xlim])
% 
xc=0.1:0.01:1.9;
cnalpha=xc.*gamma((2+xc)/2)./(2.^(1-xc)*pi.*gamma(1-xc/2));
figure;
plot(xc,cnalpha);

Ncoef=3;
I=(zxf0<zxmax)&(zxf0>zxmin)&(zyf0<zymax)&(zyf0>zymin);
zxinitial=zxf0(I);
zyinitial=zyf0(I);
x=zxf(I);
y=zyf(I);
n=length(zxinitial);
A=zeros(n,(Ncoef+1)*(Ncoef+2)/2);
A(:,1)=1;
for i=1:Ncoef
    for j=1:i+1
        A(:,i*(i+1)/2+j)=zxinitial'.^(i+1-j).*zyinitial'.^(j-1);
    end
end
Bx=x'/h;
By=y'/h;
% X=(A'*A)\(A'*Bx);
% Y=(A'*A)\(A'*By);

A2=A;
NA=length(A(1,:));
posx=1:1:NA;
for k=1:NA
    X=(A'*A)\(A'*Bx);
    I=abs(X)<0.1;
    A(:,I)=[];
    posx(I)=[];
    if isempty(X(I))
        break;
    end
end
A=A2;
posy=1:1:NA;
for k=1:NA
    Y=(A'*A)\(A'*By);
    I=abs(Y)<0.1;
    A(:,I)=[];
    posy(I)=[];
    if isempty(Y(I))
        break;
    end
end

