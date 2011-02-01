%This file reads in ANSYS data and plots forces and torques vs angle of 2
%permanent magnets
%mesh size 0.003

clc

%f=textread('torqvsangle.ansys');
f=textread('sidevalues.ansys');
%f=textread('3dvalues.ansys');


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
   
elseif f(1)==2
    
    N=f(2); %number of steps
    B=6;    %number of lines until next component
    distancex=linspace(f(3),f(4),N);   %creates a vector of distances
    Fxmx=f(5:B:end);     %Force x (maxwell stress tensor)
    Fxvw=f(6:B:end);     %Force x (virtual work)
    Fymx=f(7:B:end);
    Fyvw=f(8:B:end);
    Tzmx=f(9:B:end);    %Torque z local (maxwell stress tensor)
    Tzvw=f(10:B:end);    %Torque z local (virtual work)

    %plotting the Force in x

    figure(1)

    plot(distancex,Fxmx,'b',distancex,Fxvw,'r')
    legend('maxwell stress tensor','virtual work')
    title('Force x vs distancex')
    xlabel('Distance in m')
    ylabel('Force in Newtons')

    %plotting the Force in y

    figure(2)

    plot(distancex,Fymx,'b',distancex,Fyvw,'r')
    legend('maxwell stress tensor','virtual work')
    title('Force y vs distance x')
    xlabel('Distance in m')
    ylabel('Force in Newtons')

    %plotting the Torque around magnets own centroid

    figure(3)

    plot(distancex,Tzmx,'b',distancex,Tzvw,'r')
    legend('maxwell stress tensor','virtual work')
    title('Torque Z vs distance x')
    xlabel('Distance in m')
    ylabel('Torque in Newton metres/metre')
    
elseif f(1)==3
    
    N=f(2); %number of steps
    B=8;    %number of lines until next component
    distancex=linspace(f(3),f(4),N);   %creates a vector of distances
    Fxmx=f(5:B:end);     %Force x (maxwell stress tensor)
    Fxvw=f(6:B:end);     %Force x (virtual work)
    Fymx=f(7:B:end);
    Fyvw=f(8:B:end);
    Fzmx=f(9:B:end);
    Fzvw=f(10:B:end);
    Tzmx=f(11:B:end);    %Torque z local (maxwell stress tensor)
    Tzvw=f(12:B:end);    %Torque z local (virtual work)

    %plotting the Force in x

    figure(1)

    plot(distancex,Fxmx,'b',distancex,Fxvw,'r')
    legend('maxwell stress tensor','virtual work')
    title('Force x vs distancex')
    xlabel('Distance in m')
    ylabel('Force in Newtons')

    %plotting the Force in y

    figure(2)

    plot(distancex,Fymx,'b',distancex,Fyvw,'r')
    legend('maxwell stress tensor','virtual work')
    title('Force y vs distance x')
    xlabel('Distance in m')
    ylabel('Force in Newtons')
    
    %plotting the Force in z
    
    figure(3)

    plot(distancex,Fzmx,'b',distancex,Fzvw,'r')
    legend('maxwell stress tensor','virtual work')
    title('Force z vs distance x')
    xlabel('Distance in m')
    ylabel('Force in Newtons')

    %plotting the Torque around magnets own centroid

    figure(4)

    plot(distancex,Tzmx,'b',distancex,Tzvw,'r')
    legend('maxwell stress tensor','virtual work')
    title('Torque Z vs distance x')
    xlabel('Distance in m')
    ylabel('Torque in Newton metres')
    
    
       
end



%%

a1 = 0.5;
b1 = 0.01;
c1 = 0.004;

a2 = 0.5;
b2 = 0.01;
c2 = 0.004;

alpha = 0;
beta = 0:0.001:0.06;
gamma = 0.024;

delta = alpha;
epsilon = beta;
zeta = gamma;

br = 1.3;

[Tx Ty Tz] = Torque(a1,b1,c1,a2,b2,c2,alpha,beta,gamma,delta,epsilon,zeta,br,-br);

figure(5); hold on;

plot(beta,real(Tx),'r');
plot(beta,real(Ty),'b');
plot(beta,real(Tz),'k');
legend('x','y','z')


%% check forces


magnet_fixed.dim = [0.02 1 0.008];
magnet_float.dim = [0.02 1 0.008];

magnet_fixed.magn = 1.3;
magnet_float.magn = 1.3;

magnet_fixed.magdir = [0 0 1]; % z
magnet_float.magdir = [0 0 -1]; % -z

N = 51;
offset = repmat([0; 0; 0.024],[1 N]);
displ = linspace(0, 0.06, N);
displ_range = offset+[1; 0; 0]*displ;

f1_xyz = magnetforces(magnet_fixed,magnet_float,displ_range);

figure(6); plot(displ,f1_xyz(3,:)) %1,2,3 here is xyz but z is vertical (y)
%and x is horizontal

