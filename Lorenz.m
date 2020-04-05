clear;
clc;

zm=2;
zxmin=-zm;
zxmax=zm;
zymin=-zm;
zymax=zm;
zzmin=-1;
zzmax=3;
Nzx=400;
Nzy=400;
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
h=0.001;
alpha=1;
sigma=2;

gamma0=h;
M=stblrnd(alpha/2,1,2*(gamma0*cos(pi*alpha/4))^(2/alpha),0,1,Nz);
Normal=randn(3,Nz);
zxf=10*(z0y-z0x)*h+sigma*sqrt(M).*Normal(1,:);
zxf0=z0x;
zyf=(4*z0x-z0y-z0x.*z0z)*h+sigma*sqrt(M).*Normal(2,:);
zyf0=z0y;
zzf=(-8/3*z0z+z0x.*z0y)*h+sigma*sqrt(M).*Normal(3,:);
zzf0=z0z;

% zxf0=zxf;
% zyf0=zyf;
xf=sqrt(zxf.^2+zyf.^2+zzf.^2);

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

cnalpha=alpha0*gamma((3+alpha0)/2)/(2^(1-alpha0)*pi^1.5*gamma(1-alpha0/2));
pos=0:N;
sigmak=(a^alpha0*m.^(alpha0*pos).*nk*alpha0/(h*Nz*4*pi*cnalpha*(1-m^(-alpha0)))).^(1/alpha0);
sigma0=sum(sigmak)/(N+1);

% figure;
% plot3(zxf(1:1e6),zyf(1:1e6),zzf(1:1e6),'.');

% xc=0.1:0.01:1.9;
% cnalpha=xc.*gamma((3+xc)/2)./(2.^(1-xc)*pi^1.5.*gamma(1-xc/2));
% figure;
% plot(xc,cnalpha);

%%% drift coefficient
Ncoef=2;
I=(xf<xlim);
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

