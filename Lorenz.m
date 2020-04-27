clear;
clc;

zm=2;
zxmin=-zm;
zxmax=zm;
zymin=-zm;
zymax=zm;
zzmin=-zm;
zzmax=zm;
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

h=0.001;
alpha=1;
sigma=2;
epsilong=1;

%%% Generate data
M=stblrnd(alpha/2,1,2*(h*cos(pi*alpha/4))^(2/alpha),0,1,Nz);
Normal=randn(3,Nz);
Bh=sqrt(h)*randn(3,Nz);
zxf=10*(z0y-z0x)*h+(1+z0z).*Bh(1,:)+Bh(2,:)+sigma*sqrt(M).*Normal(1,:);
zxf0=z0x;
zyf=(4*z0x-z0y-z0x.*z0z)*h+z0y.*Bh(2,:)+sigma*sqrt(M).*Normal(2,:);
zyf0=z0y;
zzf=(-8/3*z0z+z0x.*z0y)*h+z0x.*Bh(3,:)+sigma*sqrt(M).*Normal(3,:);
zzf0=z0z;

xf=sqrt(zxf.^2+zyf.^2+zzf.^2);

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

cnalpha=alpha0*gamma((3+alpha0)/2)/(2^(1-alpha0)*pi^1.5*gamma(1-alpha0/2));
pos=0:N;
sigmak=(q^alpha0*m.^(alpha0*pos).*nk*alpha0/(h*Nz*4*pi*cnalpha*(1-m^(-alpha0)))).^(1/alpha0);
sigma0=sum(sigmak)/(N+1);

% figure;
% plot3(zxf(1:1e6),zyf(1:1e6),zzf(1:1e6),'.');

% xc=0.1:0.01:1.9;
% cnalpha=xc.*gamma((3+xc)/2)./(2.^(1-xc)*pi^1.5.*gamma(1-xc/2));
% figure;
% plot(xc,cnalpha);

% %%% Identify drift term
% Ncoef=2;
% I=(xf<epsilong);
% zxinitial=zxf0(I);
% zyinitial=zyf0(I);
% zzinitial=zzf0(I);
% x=zxf(I);
% y=zyf(I);
% z=zzf(I);
% n=length(zxinitial);
% A=zeros(n,(Ncoef+1)*(Ncoef+2)*(Ncoef+3)/6);
% A(:,1)=1;
% for i=1:Ncoef
%     for j=i:(-1):0
%         for k=i-j:(-1):0
%             A(:,i*(i+1)*(i+2)/6+(i-j)*(i-j+1)/2+i-j-k+1)=zxinitial'.^j.*zyinitial'.^k.*zzinitial'.^(i-j-k);
%         end
%     end
% end
% A2=A;
% Bx=n/Nz*x'/h;
% By=n/Nz*y'/h;
% Bz=n/Nz*z'/h;
% 
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
% 
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
% A=A2;
% posz=1:1:NA;
% for k=1:NA
%     Z=(A'*A)\(A'*Bz);
%     I=abs(Z)<0.05;
%     A(:,I)=[];
%     posz(I)=[];
%     if isempty(Z(I))
%         break;
%     end
% end
% 
% %%% Identify diffusion term
% Bx=n/Nz*x.^2'/h-4/3*pi*sigma0^alpha0*epsilong^(2-alpha0)*cnalpha/(2-alpha0);
% By=n/Nz*y.^2'/h-4/3*pi*sigma0^alpha0*epsilong^(2-alpha0)*cnalpha/(2-alpha0);
% Bz=n/Nz*z.^2'/h-4/3*pi*sigma0^alpha0*epsilong^(2-alpha0)*cnalpha/(2-alpha0);
% Bxy=n/Nz*(x.*y)'/h;
% Bxz=n/Nz*(x.*z)'/h;
% Byz=n/Nz*(y.*z)'/h;
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
% posBtz=1:1:NA;
% for k=1:NA
%     X33=(A'*A)\(A'*Bz);
%     I=abs(X33)<0.05;
%     A(:,I)=[];
%     posBtz(I)=[];
%     if isempty(X33(I))
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
% A=A2;
% posBtxz=1:1:NA;
% for k=1:NA
%     X13=(A'*A)\(A'*Bxz);
%     I=abs(X13)<0.05;
%     A(:,I)=[];
%     posBtxz(I)=[];
%     if isempty(X13(I))
%         break;
%     end
% end
% A=A2;
% posBtyz=1:1:NA;
% for k=1:NA
%     X23=(A'*A)\(A'*Byz);
%     I=abs(X23)<0.05;
%     A(:,I)=[];
%     posBtyz(I)=[];
%     if isempty(X23(I))
%         break;
%     end
% end

%%% Identify drift and diffusion terms using SSR sparsing
Ncoef=2;
I=(xf<epsilong);
zxinitial=zxf0(I);
zyinitial=zyf0(I);
zzinitial=zzf0(I);
x=zxf(I);
y=zyf(I);
z=zzf(I);
n=length(zxinitial);
Bx=n/Nz*x'/h;
By=n/Nz*y'/h;
Bz=n/Nz*z'/h;
Bxx=n/Nz*x.^2'/h-4/3*pi*sigma0^alpha0*epsilong^(2-alpha0)*cnalpha/(2-alpha0);
Byy=n/Nz*y.^2'/h-4/3*pi*sigma0^alpha0*epsilong^(2-alpha0)*cnalpha/(2-alpha0);
Bzz=n/Nz*z.^2'/h-4/3*pi*sigma0^alpha0*epsilong^(2-alpha0)*cnalpha/(2-alpha0);
Bxy=n/Nz*(x.*y)'/h;
Bxz=n/Nz*(x.*z)'/h;
Byz=n/Nz*(y.*z)'/h;
nbinsq=10;
nbin=nbinsq^3;
xbin=zeros(1,nbin);
ybin=zeros(1,nbin);
zbin=zeros(1,nbin);
Bxbin=zeros(1,nbin);
Bybin=zeros(1,nbin);
Bzbin=zeros(1,nbin);
Bbinxx=zeros(1,nbin);
Bbinyy=zeros(1,nbin);
Bbinzz=zeros(1,nbin);
Bbinxy=zeros(1,nbin);
Bbinxz=zeros(1,nbin);
Bbinyz=zeros(1,nbin);
wbin=zeros(1,nbin);
xl=linspace(zxmin-0.001,zxmax+0.001,nbinsq+1);
yl=linspace(zymin-0.001,zymax+0.001,nbinsq+1);
zl=linspace(zzmin-0.001,zzmax+0.001,nbinsq+1);
for k=1:nbinsq
    for j=1:nbinsq
        for i=1:nbinsq
            Ix=(zxinitial>xl(i))&(zxinitial<xl(i+1));
            Iy=(zyinitial>yl(j))&(zyinitial<yl(j+1));
            Iz=(zzinitial>zl(k))&(zzinitial<zl(k+1));
            I=Ix&Iy&Iz;
            n0=length(zxinitial(I));
            wbin((k-1)*nbinsq^2+(j-1)*nbinsq+i)=n0;
            xbin((k-1)*nbinsq^2+(j-1)*nbinsq+i)=sum(zxinitial(I))/n0;
            ybin((k-1)*nbinsq^2+(j-1)*nbinsq+i)=sum(zyinitial(I))/n0;
            zbin((k-1)*nbinsq^2+(j-1)*nbinsq+i)=sum(zzinitial(I))/n0;
            Bxbin((k-1)*nbinsq^2+(j-1)*nbinsq+i)=sum(Bx(I))/n0;
            Bybin((k-1)*nbinsq^2+(j-1)*nbinsq+i)=sum(By(I))/n0;
            Bzbin((k-1)*nbinsq^2+(j-1)*nbinsq+i)=sum(Bz(I))/n0;
            Bbinxx((k-1)*nbinsq^2+(j-1)*nbinsq+i)=sum(Bxx(I))/n0;
            Bbinyy((k-1)*nbinsq^2+(j-1)*nbinsq+i)=sum(Byy(I))/n0;
            Bbinzz((k-1)*nbinsq^2+(j-1)*nbinsq+i)=sum(Bzz(I))/n0;
            Bbinxy((k-1)*nbinsq^2+(j-1)*nbinsq+i)=sum(Bxy(I))/n0;
            Bbinxz((k-1)*nbinsq^2+(j-1)*nbinsq+i)=sum(Bxz(I))/n0;
            Bbinyz((k-1)*nbinsq^2+(j-1)*nbinsq+i)=sum(Byz(I))/n0;
        end
    end
end
wbin=diag(wbin)/(n/nbin);
A=zeros(nbin,(Ncoef+1)*(Ncoef+2)*(Ncoef+3)/6);
A(:,1)=1;
for i=1:Ncoef
    for j=i:(-1):0
        for k=i-j:(-1):0
            A(:,i*(i+1)*(i+2)/6+(i-j)*(i-j+1)/2+i-j-k+1)=xbin'.^j.*ybin'.^k.*zbin'.^(i-j-k);
        end
    end
end
A=wbin*A;
Bx=wbin*Bxbin';
By=wbin*Bybin';
Bz=wbin*Bzbin';
Bxx=wbin*Bbinxx';
Byy=wbin*Bbinyy';
Bzz=wbin*Bbinzz';
Bxy=wbin*Bbinxy';
Bxz=wbin*Bbinxz';
Byz=wbin*Bbinyz';

% non-sparse solution
X0=(A'*A)\(A'*Bx);
Y0=(A'*A)\(A'*By);
Z0=(A'*A)\(A'*Bz);
XX0=(A'*A)\(A'*Bxx);
YY0=(A'*A)\(A'*Byy);
ZZ0=(A'*A)\(A'*Bzz);
XY0=(A'*A)\(A'*Bxy);
XZ0=(A'*A)\(A'*Bxz);
YZ0=(A'*A)\(A'*Byz);

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
delta0=3;
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

% drift b3(x)
A2=A3;
Bz2=Bz;
deleteorderz=zeros(1,NA);
deltaSSRz=zeros(1,NA);
posz=1:1:NA;
for k=1:NA
    A=A2;
    Bz=Bz2;
    Z=(A'*A)\(A'*Bz);
    [m,I]=min(abs(Z));
    deleteorderz(k)=posz(I);
    A(:,I)=[];
    posz(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRz(k)=norm(Bz)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        Bz=Bz2;
        Ai=A(IA(i,:)',:);
        Bi=Bz(IA(i,:)');
        A(IA(i,:)',:)=[];
        Bz(IA(i,:)')=[];
        SX=(A'*A)\(A'*Bz);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRz(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRz(k)/deltaSSRz(k-1)>delta0)&&(deltaSSRz(k-1)/deltaSSRz(k-2)<delta0)
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

% diffusion sigma33(x)
A2=A3;
Bzz2=Bzz;
deleteorderzz=zeros(1,NA);
deltaSSRzz=zeros(1,NA);
poszz=1:1:NA;
for k=1:NA
    A=A2;
    Bzz=Bzz2;
    ZZ=(A'*A)\(A'*Bzz);
    [m,I]=min(abs(ZZ));
    deleteorderzz(k)=poszz(I);
    A(:,I)=[];
    poszz(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRzz(k)=norm(Bzz)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        Bzz=Bzz2;
        Ai=A(IA(i,:)',:);
        Bi=Bzz(IA(i,:)');
        A(IA(i,:)',:)=[];
        Bzz(IA(i,:)')=[];
        SX=(A'*A)\(A'*Bzz);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRzz(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRzz(k)/deltaSSRzz(k-1)>delta0)&&(deltaSSRzz(k-1)/deltaSSRzz(k-2)<delta0)
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

% diffusion sigma13(x)
A2=A3;
Bxz2=Bxz;
deleteorderxz=zeros(1,NA);
deltaSSRxz=zeros(1,NA);
posxz=1:1:NA;
for k=1:NA
    A=A2;
    Bxz=Bxz2;
    XZ=(A'*A)\(A'*Bxz);
    [m,I]=min(abs(XZ));
    deleteorderxz(k)=posxz(I);
    A(:,I)=[];
    posxz(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRxz(k)=norm(Bxz)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        Bxz=Bxz2;
        Ai=A(IA(i,:)',:);
        Bi=Bxz(IA(i,:)');
        A(IA(i,:)',:)=[];
        Bxz(IA(i,:)')=[];
        SX=(A'*A)\(A'*Bxz);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRxz(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRxz(k)/deltaSSRxz(k-1)>delta0)&&(deltaSSRxz(k-1)/deltaSSRxz(k-2)<delta0)
        break;
    end
end

% diffusion sigma23(x)
A2=A3;
Byz2=Byz;
deleteorderyz=zeros(1,NA);
deltaSSRyz=zeros(1,NA);
posyz=1:1:NA;
for k=1:NA
    A=A2;
    Byz=Byz2;
    YZ=(A'*A)\(A'*Byz);
    [m,I]=min(abs(YZ));
    deleteorderyz(k)=posyz(I);
    A(:,I)=[];
    posyz(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRyz(k)=norm(Byz)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        Byz=Byz2;
        Ai=A(IA(i,:)',:);
        Bi=Byz(IA(i,:)');
        A(IA(i,:)',:)=[];
        Byz(IA(i,:)')=[];
        SX=(A'*A)\(A'*Byz);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRyz(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRyz(k)/deltaSSRyz(k-1)>delta0)&&(deltaSSRyz(k-1)/deltaSSRyz(k-2)<delta0)
        break;
    end
end


