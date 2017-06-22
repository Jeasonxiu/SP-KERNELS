clear;

the=10*pi/180;
beta=4.4;
alpha=7.920;

%beta=1.0;
%alpha=1.0;

p=sin(the)/beta;

qa0=sqrt(alpha^-2 - p^2);
qb0=sqrt(beta ^-2 - p^2);

Vpz=-(1.-2.*beta^2*p^2) / (2.*alpha*qa0);
Vpr=p*beta^2 / alpha;
Vsz=p*beta;
Vsr=(1.-2.*beta^2*p^2) / (2.*beta*qb0);

M=zeros(2,2);
M(1,1)=Vpz;
M(1,2)=Vpr;
M(2,1)=Vsz;
M(2,2)=Vsr;

[V,D]=eig(M);
disp(V)

figure(1); clf;
quiver(0,0,V(1,1),V(2,1),'r'); hold on;
quiver(0,0,V(1,2),V(2,2),'b')
ylim([-1,1]); xlim([-1,1])