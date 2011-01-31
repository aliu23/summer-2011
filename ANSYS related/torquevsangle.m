%This file reads in ANSYS data and plots forces and torques vs angle of 2
%permanent magnets
%mesh size 0.003

clc

f=textread('torqvsangle.ansys');

if f(1)==1
    
   N=f(2); %number of steps
   B=8;    %number of lines until next component
   angles=linspace((f(3).*180/pi),f(4).*180/pi,N);   %creates a vector of angles 
   Fxmx=f(5:B:end);     %Force x (maxwell stress tensor)
   Fxvw=f(6:B:end);     %Force x (virtual work)
   Fymx=f(7:B:end);
   Fyvw=f(8:B:end);
   Tzmx=f(11:B:end);    %Torque z local (maxwell stress tensor)
   Tzvw=f(12:B:end);    %Torque z local (virtual work)
    
    
end

%plotting the Force in x

figure(1)

plot(angles,Fxmx,'b',angles,Fxvw,'r')
legend('maxwell stress tensor','virtual work')
title('Force x vs angle')
xlabel('Angle in degrees')
ylabel('Force in Newtons')

%plotting the Force in y

figure(2)

plot(angles,Fymx,'b',angles,Fyvw,'r')
legend('maxwell stress tensor','virtual work')
title('Force y vs angle')
xlabel('Angle in degrees')
ylabel('Force in Newtons')

%plotting the Torque around magnets own centroid 

figure(3)

plot(angles,Tzmx,'b',angles,Tzvw,'r')
legend('maxwell stress tensor','virtual work')
title('Torque Z vs angle')
xlabel('Angle in degrees')
ylabel('Torque in Newton metres')




