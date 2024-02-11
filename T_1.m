clear all;
close all;
clc;

%%
U=1;
gamma=input("input the value of circulation or gamma = ");
k=input("input the doublet strength = ");
ro=1.21;
alpha=0;
a=sqrt(k/(2*pi*U));

%%
syms r theta 
z=r.*exp(theta*i)*exp(-alpha*i);
F=U*(z+a^2./z)-log(z).*(gamma/(2*pi))*i;
% W=diff(F);
W=U*(1-a^2./z.^2)-(gamma./(2*pi*z))*i;
% I=int(0.5*ro*i*W.^2,theta,[0,2*pi]);
% A=real(double(I));
% N=-imag(double(I));
B1=gamma*i*exp(-alpha*i);
G=-U*ro*B1;
A=real(G);
N=-imag(G);

B=[cos(alpha),-sin(alpha);sin(alpha),cos(alpha)]*[N;A];
L=B(1);
D=B(2);

phi=real(F);
shi=imag(F);
u_r=real(W./exp(-theta*i));
u_theta=-imag(W./exp(-theta*i));
u_theta_f=matlabFunction(u_theta);
u_r_f=matlabFunction(u_r);

%%
%location of stagnation
eqns=u_r==0;
s1=solve(eqns,[r,theta]);
r_s=double(s1.r);
eqns=u_theta_f(abs(r_s(1)),theta)==0;
s2=solve(eqns,theta);
theta_s=abs(double(s2))*180/pi;
beta=asin(gamma/(4*pi*a*U))*(180/pi);

%%
%streamline,velocity,pressure and force
shi_f=matlabFunction(shi);
p=linspace(0,3);
q=linspace(0,2*pi);
[r,theta]=meshgrid(p,q);
figure(1);
contour(r.*cos(theta),r.*sin(theta),shi_f(r,theta),101);
title("streamline");

figure(2);
t=0:pi/120:2*pi;
cp=1-4*(sin(t)-sin(beta*pi/180)).^2;
plot(t,cp);
title("pressure variation");
xlabel("--- theta --->");
ylabel("---- cp ---->");


syms t   
cp=1-4*(sin(t)-sin(beta*pi/180)).^2;
J1=matlabFunction(cp.*sin(t));
J2=matlabFunction(cp.*cos(t));
I1=int(J1,t,[0,2*pi]);
I2=int(J2,t,[0,2*pi]);
Fy=0.5*ro*a*U^2*double(I1);
Fx=0.5*ro*a*U^2*double(I2);

figure(3);
t=0:pi/120:2*pi;
v=sqrt((u_r_f(a,t)).^2+(u_theta_f(a,t)).^2);
plot(t,v);
title("velocity");
xlabel("----- theta --->");
ylabel("----- velocity --->");

