clear;
clc;

zmin=-3;
zmax=3;
Nz=1e7;
z0=linspace(zmin,zmax,Nz);

% N=1;                   %%% 一个初始点给出的样本数
h=0.001;
alpha=0.5;
sigma=2;
epsilong=1;

%%% Generate data
M=h^(1/alpha)*stblrnd(alpha,0,1,0,1,Nz);
Bh=sqrt(h)*randn(1,Nz);
xf=(4*z0-z0.^3)*h+(1+z0).*Bh+sigma*M;
xf0=z0;

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
% nratio=nk(1:end-1)./nk(2:end);
% alpha1=log(nratio)./(log(m));
% alpha0=sum(alpha1)/N;

cnalpha=alpha0*gamma((1+alpha0)/2)/(2^(1-alpha0)*sqrt(pi)*gamma(1-alpha0/2));
pos=0:N;
sigmak=q*m.^pos.*(nk*alpha0/(h*Nz*2*cnalpha*(1-m^(-alpha0)))).^(1/alpha0);
sigma0=sum(sigmak)/(N+1);

% xc=0.1:0.01:1.9;
% cnalpha=xc.*gamma((1+xc)/2)./(2.^(1-xc)*sqrt(pi).*gamma(1-xc/2));
% figure;
% plot(xc,cnalpha);
% 
% %%% Identify drift term
% Ncoef=6;
% I=(abs(xf)<epsilong);
% zinitial=xf0(I);
% x=xf(I);
% n=length(zinitial);
% A=zeros(n,Ncoef+1);
% A(:,1)=1;
% for i=1:Ncoef
%     A(:,i+1)=zinitial'.^i;
% end
% A2=A;
% B=n/Nz*x'/h;
% % X=(A'*A)\(A'*B);
% pos=0:1:Ncoef;
% for k=1:Ncoef
%     X=(A'*A)\(A'*B);
%     I=abs(X)<0.05;
%     A(:,I)=[];
%     pos(I)=[];
%     if isempty(X(I))
%         break;
%     end
% end
% 
% %%% Identify diffusion term
% A=A2;
% B=n/Nz*x.^2'/h-2*sigma0^alpha0*epsilong^(2-alpha0)*cnalpha/(2-alpha0);
% % XBt=(A'*A)\(A'*B);
% posBt=0:1:Ncoef;
% for k=1:Ncoef
%     XBt=(A'*A)\(A'*B);
%     I=abs(XBt)<0.05;
%     A(:,I)=[];
%     posBt(I)=[];
%     if isempty(XBt(I))
%         break;
%     end
% end

%%% Identify drift and diffusion terms using SSR sparsing
Ncoef=6;
I=(abs(xf)<epsilong);
zinitial=xf0(I);
x=xf(I);
n=length(zinitial);
B=n/Nz*x'/h;
BBt=n/Nz*x.^2'/h-2*sigma0^alpha0*epsilong^(2-alpha0)*cnalpha/(2-alpha0);
nbin=1e3;
xbin=zeros(1,nbin);
Bbin=zeros(1,nbin);
BbinBt=zeros(1,nbin);
wbin=zeros(1,nbin);
n0=floor(n/nbin);
wbin(1:end-1)=n0;
wbin(end)=(n-n0*(nbin-1));
xbin(1)=sum(zinitial(1:n0))/n0;
Bbin(1)=sum(B(1:n0))/n0;
BbinBt(1)=sum(B(1:n0))/n0;
for i=2:nbin
    I=sum(wbin(1:i-1))+1:sum(wbin(1:i));
    xbin(i)=sum(zinitial(I))/length(I);
    Bbin(i)=sum(B(I))/length(I);
    BbinBt(i)=sum(BBt(I))/length(I);
end
wbin=diag(wbin)/n0;
A=zeros(nbin,Ncoef+1);
A(:,1)=1;
for i=1:Ncoef
    A(:,i+1)=xbin'.^i;
end
A0=wbin*A;
B0=wbin*Bbin';
B0Bt=wbin*BbinBt';

X0=(A0'*A0)\(A0'*B0);
XBt0=(A0'*A0)\(A0'*B0Bt);

kfold=2;
nblock=floor(nbin/kfold);
n=kfold*nblock;
A=A0(1:nbin,:);
B=B0(1:nbin);
BBt=B0Bt(1:nbin);
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

A3=A;
A2=A;
B2=B;
deleteorder=zeros(1,Ncoef+1);
deltaSSR=zeros(1,Ncoef+1);
delta0=1.5;
pos=0:1:Ncoef;
for k=1:Ncoef+1
    A=A2;
    B=B2;
    X=(A'*A)\(A'*B);
    [m,I]=min(abs(X));
    deleteorder(k)=pos(I);
    A(:,I)=[];
    pos(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSR(k)=norm(B)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        B=B2;
        Ai=A(IA(i,:)',:);
        Bi=B(IA(i,:)');
        A(IA(i,:)',:)=[];
        B(IA(i,:)')=[];
        SX=(A'*A)\(A'*B);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSR(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSR(k)/deltaSSR(k-1)>delta0)&&(deltaSSR(k-1)/deltaSSR(k-2)<delta0)
        break;
    end
end

A2=A3;
B2=BBt;
deleteorderBt=zeros(1,Ncoef+1);
deltaSSRBt=zeros(1,Ncoef+1);
posBt=0:1:Ncoef;
for k=1:Ncoef+1
    A=A2;
    B=B2;
    XBt=(A'*A)\(A'*B);
    [m,I]=min(abs(XBt));
    deleteorderBt(k)=posBt(I);
    A(:,I)=[];
    posBt(I)=[];
    A2=A;
    
    if isempty(A)
        deltaSSRBt(k)=norm(B)/kfold;
        break;
    end
    
    SSR=0;
    for i=1:kfold
        A=A2;
        B=B2;
        Ai=A(IA(i,:)',:);
        Bi=B(IA(i,:)');
        A(IA(i,:)',:)=[];
        B(IA(i,:)')=[];
        SX=(A'*A)\(A'*B);
        SSR=SSR+norm(Bi-Ai*SX)^2;
    end
    deltaSSRBt(k)=sqrt(SSR/kfold);
    
    if (k>2)&&(deltaSSRBt(k)/deltaSSRBt(k-1)>delta0)&&(deltaSSRBt(k-1)/deltaSSRBt(k-2)<delta0)
        break;
    end
end
