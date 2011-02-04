%% 3D analysis of sliding parallel magnets
%
% This file reads in ANSYS data for 2 parallel magnets and plots 
% Forces and Torques vs a horizontal displacement of the top magnet

%% Setup

clc

data=readansys('3dvalues.ansys');

%% 

figure(7); hold on

plot(data.distance,data.force_mx(:,1),'b')
plot(data.distance,data.force_vw(:,1),'r')
legend('maxwell stress tensor','virtual work')
title('Force  x vs Distance x')
xlabel('Distance x (m)')
ylabel('Force (N)')


figure(8); hold on

plot(data.distance,data.force_mx(:,2),'b')
plot(data.distance,data.force_vw(:,2),'r')
legend('maxwell stress tensor','virtual work')
title('Force y vs Distance x')
xlabel('Distance x (m)')
ylabel('Force (N)')

figure(9); hold on

plot(data.distance,data.force_mx(:,3),'b')
plot(data.distance,data.force_vw(:,3),'r')
legend('maxwell stress tensor','virtual work')
title('Force z vs Distance x')
xlabel('Distance x (m)')
ylabel('Force (N)')


figure(10); hold on

plot(data.distance,data.torque_mx,'b')
plot(data.distance,data.torque_vw,'r')
legend('maxwell stress tensor','virtual work')
title('Torque z vs distance x')
xlabel('Distance x (m)')
ylabel('Torque (N*m)')

