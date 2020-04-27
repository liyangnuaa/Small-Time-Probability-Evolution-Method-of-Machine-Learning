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

% N=1;                   %%% 一个初始点给出的样本数
h=0.001;
alpha=1;
sigma=2;
epsilong=1;

%%% Generate data
M=stblrnd(alpha/2,1,2*(h*cos(pi*alpha/4))^(2/alpha),0,1,Nz);
Normal=randn(2,Nz);
Bh=sqrt(h)*randn(2,Nz);
zxf=(z0x-z0x.^3-5*z0x.*z0y.^2)*h+(1+z0y).*Bh(1,:)+Bh(2,:)+sigma*sqrt(M).*Normal(1,:);
zxf0=z0x;
zyf=-(1+z0x.^2).*z0y*h+z0x.*Bh(2,:)+sigma*sqrt(M).*Normal(2,:);
zyf0=z0y;
clear z0x z0y

xf=sqrt(zxf.^2+zyf.^2);

%%% Identify alpha and sigma
q=epsilong;
m=5;
N=2;
nk=zeros(1,N+1);
for k=0:N
    I=(abs(xf)>=m^k*q)&(abs(xf)<m^(k+1)*q);
    nk(k+1)=length(xf(I));
end
nratio=nk(1)./nk(2:end);
pos=1:N;
alpha1=log(nratio)./(pos*log(m));
alpha0=sum(alpha1)/N;

cnalpha=alpha0*gamma((2+alpha0)/2)/(2^(1-alpha0)*pi*gamma(1-alpha0/2));
pos=0:N;
sigmak=(q^alpha0*m.^(alpha0*pos).*nk*alpha0/(h*Nz*2*pi*cnalpha*(1-m^(-alpha0)))).^(1/alpha0);
sigma0=sum(sigmak)/(N+1);

% figure;
% plot(zxf0(1:1e7),zyf0(1:1e7),'.');
% axis([-xlim xlim -xlim xlim])
% 
% xc=0.1:0.01:1.9;
% cnalpha=xc.*gamma((2+xc)/2)./(2.^(1-xc)*pi.*gamma(1-xc/2));
% figure;
% plot(xc,cnalpha);

% %%% Identify drift term
% Ncoef=3;
% I=(xf<epsilong);
% zxinitial=zxf0(I);
% zyinitial=zyf0(I);
% x=zxf(I);
% y=zyf(I);
% n=length(zxinitial);
% A=zeros(n,(Ncoef+1)*(Ncoef+2)/2);
% A(:,1)=1;
% for i=1:Ncoef
%     for j=1:i+1
%         A(:,i*(i+1)/2+j)=zxinitial'.^(i+1-j).*zyinitial'.^(j-1);
%     end
% end
% Bx=n/Nz*x'/h;
% By=n/Nz*y'/h;
% % X=(A'*A)\(A'*Bx);
% % Y=(A'*A)\(A'*By);
% 
% A2=A;
% NA=length(A(1,:));
% posx=1:1:NA;
% for k=1:NA
%     X=(A'*A)\(A'*Bx);
%     I=abs(X)<0.05;
%     A(:,I)=[];
%     posx(I)=[];
%     if isempty(X(I))
%         break;
%     end
% end
% A=A2;
% posy=1:1:NA;
% for k=1:NA
%     Y=(A'*A)\(A'*By);
%     I=abs(Y)<0.05;
%     A(:,I)=[];
%     posy(I)=[];
%     if isempty(Y(I))
%         break;
%     end
% end
% 
% %%% Identify diffusion term
% Bx=n/Nz*x.^2'/h-pi*sigma0^alpha0*epsilong^(2-alpha0)*cnalpha/(2-alpha0);
% By=n/Nz*y.^2'/h-pi*sigma0^alpha0*epsilong^(2-alpha0)*cnalpha/(2-alpha0);
% Bxy=n/Nz*(x.*y)'/h;
% % X11=(A'*A)\(A'*Bx);
% % X22=(A'*A)\(A'*By);
% % X12=(A'*A)\(A'*Bxy);
% 
% A=A2;
% posBtx=1:1:NA;
% for k=1:NA
%     X11=(A'*A)\(A'*Bx);
%     I=abs(X11)<0.05;
%     A(:,I)=[];
%     posBtx(I)=[];
%     if isempty(X11(I))
%         break;
%     end
% end
% A=A2;
% posBty=1:1:NA;
% for k=1:NA
%     X22=(A'*A)\(A'*By);
%     I=abs(X22)<0.05;
%     A(:,I)=[];
%     posBty(I)=[];
%     if isempty(X22(I))
%         break;
%     end
% end
% A=A2;
% posBtxy=1:1:NA;
% for k=1:NA
%     X12=(A'*A)\(A'*Bxy);
%     I=abs(X12)<0.05;
%     A(:,I)=[];
%     posBtxy(I)=[];
%     if isempty(X12(I))
%         break;
%     end
% end

%%% Identify drift and diffusion terms using SSR sparsing
Ncoef=3;
I=(xf<epsilong);
zxinitial=zxf0(I);
zyinitial=zyf0(I);
x=zxf(I);
y=zyf(I);
n=length(zxinitial);
Bx=n/Nz*x'/h;
By=n/Nz*y'/h;
Bxx=n/Nz*x.^2'/h-pi*sigma0^alpha0*epsilong^(2-alpha0)*cnalpha/(2-alpha0);
Byy=n/Nz*y.^2'/h-pi*sigma0^alpha0*epsilong^(2-alpha0)*cnalpha/(2-alpha0);
Bxy=n/Nz*(x.*y)'/h;
nbinsq=30;
nbin=nbinsq^2;
xbin=zeros(1,nbin);
ybin=zeros(1,nbin);
Bxbin=zeros(1,nbin);
Bybin=zeros(1,nbin);
Bbinxx=zeros(1,nbin);
Bbinyy=zeros(1,nbin);
Bbinxy=zeros(1,nbin);
wbin=zeros(1,nbin);
xl=linspace(zxmin-0.001,zxmax+0.001,nbinsq+1);
yl=linspace(zymin-0.001,zymax+0.001,nbinsq+1);
for j=1:nbinsq
    for i=1:nbinsq
        Ix=(zxinitial>xl(i))&(zxinitial<xl(i+1));
        Iy=(zyinitial>yl(j))&(zyinitial<yl(j+1));
        I=Ix&Iy;
        n0=length(zxinitial(I));
        wbin((j-1)*nbinsq+i)=n0;
        xbin((j-1)*nbinsq+i)=sum(zxinitial(I))/n0;
        ybin((j-1)*nbinsq+i)=sum(zyinitial(I))/n0;
        Bxbin((j-1)*nbinsq+i)=sum(Bx(I))/n0;
        Bybin((j-1)*nbinsq+i)=sum(By(I))/n0;
        Bbinxx((j-1)*nbinsq+i)=sum(Bxx(I))/n0;
        Bbinyy((j-1)*nbinsq+i)=sum(Byy(I))/n0;
        Bbinxy((j-1)*nbinsq+i)=sum(Bxy(I))/n0;
    end
end
wbin=diag(wbin)/(n/nbin);
A=zeros(nbin,(Ncoef+1)*(Ncoef+2)/2);
A(:,1)=1;
for i=1:Ncoef
    for j=1:i+1
        A(:,i*(i+1)/2+j)=xbin'.^(i+1-j).*ybin'.^(j-1);
    end
end
A=wbin*A;
Bx=wbin*Bxbin';
By=wbin*Bybin';
Bxx=wbin*Bbinxx';
Byy=wbin*Bbinyy';
Bxy=wbin*Bbinxy';

X0=(A'*A)\(A'*Bx);
Y0=(A'*A)\(A'*By);
XX0=(A'*A)\(A'*Bxx);
YY0=(A'*A)\(A'*Byy);
XY0=(A'*A)\(A'*Bxy);

kfold=2;
nblock=floor(nbin/kfold);
n=kfold*nblock;
IA=zeros(kfold,nblock);
position=1:n;
for i=1:kfold-1
    u=rand(1,nblock);
    J=kfold+1-i;
    p=0:J:J*(nblock-1);
    IA(i,:)=position(p+floor(J*u)+1);
    position(p+floor(J*u)+1)=[];
end
IA(kfold,:)=position;

% drift b1(x)
A3=A;
A2=A;
Bx2=Bx;
NA=length(A(1,:));
deleteorderx=zeros(1,NA);
deltaSSRx=zeros(1,NA);
delta0=2;
posx=1:1:NA;
for k=1:NA
    A=A2;
    Bx=Bx2;
    X=(A'*A)\(A'*Bx);
    [m,I]=min(abs(X));
    deleteorderx(k)=posx(I);
    A(:,I)=[];
    posx(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRx(k)=norm(Bx)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        Bx=Bx2;
        Ai=A(IA(i,:)',:);
        Bi=Bx(IA(i,:)');
        A(IA(i,:)',:)=[];
        Bx(IA(i,:)')=[];
        SX=(A'*A)\(A'*Bx);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRx(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRx(k)/deltaSSRx(k-1)>delta0)&&(deltaSSRx(k-1)/deltaSSRx(k-2)<delta0)
        break;
    end
end

% drift b2(x)
A2=A3;
By2=By;
deleteordery=zeros(1,NA);
deltaSSRy=zeros(1,NA);
posy=1:1:NA;
for k=1:NA
    A=A2;
    By=By2;
    Y=(A'*A)\(A'*By);
    [m,I]=min(abs(Y));
    deleteordery(k)=posy(I);
    A(:,I)=[];
    posy(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRy(k)=norm(By)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        By=By2;
        Ai=A(IA(i,:)',:);
        Bi=By(IA(i,:)');
        A(IA(i,:)',:)=[];
        By(IA(i,:)')=[];
        SX=(A'*A)\(A'*By);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRy(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRy(k)/deltaSSRy(k-1)>delta0)&&(deltaSSRy(k-1)/deltaSSRy(k-2)<delta0)
        break;
    end
end

% diffusion sigma11(x)
A2=A3;
Bxx2=Bxx;
deleteorderxx=zeros(1,NA);
deltaSSRxx=zeros(1,NA);
posxx=1:1:NA;
for k=1:NA
    A=A2;
    Bxx=Bxx2;
    XX=(A'*A)\(A'*Bxx);
    [m,I]=min(abs(XX));
    deleteorderxx(k)=posxx(I);
    A(:,I)=[];
    posxx(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRxx(k)=norm(Bxx)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        Bxx=Bxx2;
        Ai=A(IA(i,:)',:);
        Bi=Bxx(IA(i,:)');
        A(IA(i,:)',:)=[];
        Bxx(IA(i,:)')=[];
        SX=(A'*A)\(A'*Bxx);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRxx(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRxx(k)/deltaSSRxx(k-1)>delta0)&&(deltaSSRxx(k-1)/deltaSSRxx(k-2)<delta0)
        break;
    end
end

% diffusion sigma22(x)
A2=A3;
Byy2=Byy;
deleteorderyy=zeros(1,NA);
deltaSSRyy=zeros(1,NA);
posyy=1:1:NA;
for k=1:NA
    A=A2;
    Byy=Byy2;
    YY=(A'*A)\(A'*Byy);
    [m,I]=min(abs(YY));
    deleteorderyy(k)=posyy(I);
    A(:,I)=[];
    posyy(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRyy(k)=norm(Byy)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        Byy=Byy2;
        Ai=A(IA(i,:)',:);
        Bi=Byy(IA(i,:)');
        A(IA(i,:)',:)=[];
        Byy(IA(i,:)')=[];
        SX=(A'*A)\(A'*Byy);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRyy(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRyy(k)/deltaSSRyy(k-1)>delta0)&&(deltaSSRyy(k-1)/deltaSSRyy(k-2)<delta0)
        break;
    end
end

% diffusion sigma12(x)
A2=A3;
Bxy2=Bxy;
deleteorderxy=zeros(1,NA);
deltaSSRxy=zeros(1,NA);
posxy=1:1:NA;
for k=1:NA
    A=A2;
    Bxy=Bxy2;
    XY=(A'*A)\(A'*Bxy);
    [m,I]=min(abs(XY));
    deleteorderxy(k)=posxy(I);
    A(:,I)=[];
    posxy(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRxy(k)=norm(Bxy)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        Bxy=Bxy2;
        Ai=A(IA(i,:)',:);
        Bi=Bxy(IA(i,:)');
        A(IA(i,:)',:)=[];
        Bxy(IA(i,:)')=[];
        SX=(A'*A)\(A'*Bxy);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRxy(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRxy(k)/deltaSSRxy(k-1)>delta0)&&(deltaSSRxy(k-1)/deltaSSRxy(k-2)<delta0)
        break;
    end
end
