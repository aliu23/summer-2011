
%% check torques

%This section uses a user written matlab function Torque.m in order to
%evaluate the torques between 2 parallel magnets

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

figure(11); hold on;

plot(beta,real(Tx),'k--'); %Tx is actually Tz refer to coordinates
plot(data.distance,data.torque_mx,'b')
plot(data.distance,data.torque_vw,'r')
title('Torque vs distance')
xlabel('distance (m)')
ylabel('Torque (N*m)')
legend('matlab function','ANSYS (Maxwell)','ANSYS (Virtual Work)')


%% check forces

%This section uses a user written matlab function magnetforces.m in order
%to evaluate the forces in the x y and z directions between 2 parallel
%magnets


magnet_fixed.dim = [0.02 1 0.008];
magnet_float.dim = [0.02 1 0.008];

magnet_fixed.magn = 1.3;
magnet_float.magn = 1.3;

magnet_fixed.magdir = [0 0 1]; % z
magnet_float.magdir = [0 0 -1]; % -z

N = 51;
offset = repmat([0; 0; 0.024],[1 N]); %distance between magnet centres
displ = linspace(0, 0.06, N); % distance vector
displ_range = offset+[1; 0; 0]*displ;

f1_xyz = magnetforces(magnet_fixed,magnet_float,displ_range);

figure(12); hold on 
plot(displ,f1_xyz(1,:),'k--') %1,2,3 here is xyz but z is vertical (y)
%and x is horizontal
plot(data.distance,data.force_mx(:,1),'b')
plot(data.distance,data.force_vw(:,1),'r')
title('Force x vs distance x')
xlabel('distance x (m)')
ylabel('Force (N)')
legend('matlab function','maxwell stress tensor','virtual work')

figure(13); hold on
plot(displ,f1_xyz(3,:),'k--')
plot(data.distance,data.force_mx(:,2),'b')
plot(data.distance,data.force_vw(:,2),'r')
title('Force y vs distance x')
xlabel('distance x (m)')
ylabel('Force (N)')
legend('matlab function','maxwell stress tensor','virtual work')

figure(14); hold on
plot(displ,f1_xyz(2,:))
title('Force z vs distance x')
xlabel('distance x (m)')
ylabel('Force (N)')
%legend('matlab function','maxwell stress tensor','virtual work')


