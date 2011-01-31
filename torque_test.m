%%

clear all
close all
clc

%%

a1 = 0.005;
b1 = 0.013;
c1 = 0.007;

a2 = 0.007;
b2 = 0.013;
c2 = 0.005;

alpha = 0:0.001:0.035;
beta = -0.008;
gamma = 0.015;

delta = 0;
epsilon = 0;
zeta = -0.047;

br = 1.23;

[Tx Ty Tz] = Torque(a1,b1,c1,a2,b2,c2,alpha,beta,gamma,delta,epsilon,zeta,br,-br);

willfig; hold on;

plot(real(Tx),'r');
plot(real(Ty),'b');
plot(real(Tz),'k');
legend('x','y','z')
