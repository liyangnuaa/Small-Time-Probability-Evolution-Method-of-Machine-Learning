clear;
clc;

zm=2;
zxmin=-zm;
zxmax=zm;
zymin=-zm;
zymax=zm;
zzmin=-1;
zzmax=3;
Nzx=500;
Nzy=500;
Nzz=400;
Nz=Nzx*Nzy*Nzz;
y0=linspace(zymin,zymax,Nzy);
z0=linspace(zzmin,zzmax,Nzz);
z0x=zeros(1,Nz);
z0y=zeros(1,Nz);
z0z=zeros(1,Nz);
for k=1:Nzz
    for i=1:Nzy
        z0x((k-1)*Nzx*Nzy+(i-1)*Nzx+1:(k-1)*Nzx*Nzy+i*Nzx)=linspace(zxmin,zxmax,Nzx);
        z0y((k-1)*Nzx*Nzy+(i-1)*Nzx+1:(k-1)*Nzx*Nzy+i*Nzx)=y0(i);
    end
    z0z((k-1)*Nzx*Nzy+1:k*Nzx*Nzy)=z0(k);
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
zzf=zeros(1,Nz*N);
zxf0=zeros(1,Nz*N);
zyf0=zeros(1,Nz*N);
zzf0=zeros(1,Nz*N);
position=(0:Nz-1)*N;

gamma0=gamma((1+alpha)/2)*gamma(1.5)*h/(sqrt(pi)*gamma(1.5+alpha/2));
M=stblrnd(alpha/2,1,2*(gamma0*cos(pi*alpha/4))^(2/alpha),0,1,Nz+1e6);
I=abs(M)<xlim^2;
M=M(I);
M=M(1:Nz);
Normal=randn(3,Nz);
zxf(position+1)=10*(z0y-z0x)*h+sigma*sqrt(M).*Normal(1,:);
zxf0(position+1)=z0x;
zyf(position+1)=(4*z0x-z0y-z0x.*z0z)*h+sigma*sqrt(M).*Normal(2,:);
zyf0(position+1)=z0y;
zzf(position+1)=(-8/3*z0z+z0x.*z0y)*h+sigma*sqrt(M).*Normal(3,:);
zzf0(position+1)=z0z;

% zxf0=zxf;
% zyf0=zyf;
znorm=sqrt(zxf.^2+zyf.^2+zzf.^2);
xf=znorm;

xpts=linspace(1e-2,xlim,1e3);
[pdf,x00]=ksdensity(xf,xpts);
% pdf1=2*pdf/(h*2*pi);
pdf1=pdf/(h*4*pi);
x2=x00;
pdf2=pdf1;
[~,I]=min(abs(x00));
% I=1:2;
x2(I)=[];
pdf2(I)=[];
% [~,I]=min(abs(x00));
% x2(I)=[];
% pdf2(I)=[];
% [~,I]=min(abs(x00));
% x2(I)=[];
% pdf2(I)=[];
% [~,I]=min(abs(x00));
% x2(I)=[];
% pdf2(I)=[];


syms t
% f=fittype('k1/(t+r1)^r2+k2','independent','t','coefficients',{'r1','r2','k1','k2'});
% cfun=fit(x0,y0,f,'StartPoint',[3,4,-1e5,38]); %显示拟合函数，数据必须为列向量形式
f=fittype('calpha/abs(t).^(1+alpha0)','independent','t','coefficients',{'calpha','alpha0'});
cfun=fit(x2',pdf2',f,'StartPoint',[0.1,1]); %显示拟合函数，数据必须为列向量形式
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
% plot3(zxf(1:1e6),zyf(1:1e6),zzf(1:1e6),'.');

% xc=0.1:0.01:1.9;
% cnalpha=xc.*gamma((3+xc)/2)./(2.^(1-xc)*pi^1.5.*gamma(1-xc/2));
% figure;
% plot(xc,cnalpha);

%%% drift coefficient
Ncoef=2;
I=(zxf0<zxmax)&(zxf0>zxmin)&(zyf0<zymax)&(zyf0>zymin)&(zzf0<zzmax)&(zzf0>zzmin);
zxinitial=zxf0(I);
zyinitial=zyf0(I);
zzinitial=zzf0(I);
x=zxf(I);
y=zyf(I);
z=zzf(I);
n=length(zxinitial);
A=zeros(n,(Ncoef+1)*(Ncoef+2)*(Ncoef+3)/6);
A(:,1)=1;
for i=1:Ncoef
    for j=i:(-1):0
        for k=i-j:(-1):0
            A(:,i*(i+1)*(i+2)/6+(i-j)*(i-j+1)/2+i-j-k+1)=zxinitial'.^j.*zyinitial'.^k.*zzinitial'.^(i-j-k);
        end
    end
end
A2=A;
Bx=x'/h;
By=y'/h;
Bz=z'/h;

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

A=A2;
posz=1:1:NA;
for k=1:NA
    Z=(A'*A)\(A'*Bz);
    I=abs(Z)<0.1;
    A(:,I)=[];
    posz(I)=[];
    if isempty(Z(I))
        break;
    end
end

