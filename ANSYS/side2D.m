%% 2D analysis of sliding parallel magnets
%
% This file reads in ANSYS data for 2 parallel magnets and plots 
% Forces and Torques vs a horizontal displacement of the top magnet

%% Setup

clc

data=readansys('sidevalues.ansys');

%%

figure(4); hold on

plot(data.distance,data.force_mx(:,1),'b')
plot(data.distance,data.force_vw(:,1),'r')
legend('maxwell stress tensor','virtual work')
title('Force per unit length x vs Distance x')
xlabel('Distance x (m)')
ylabel('Force per unit length (N/m)')


figure(5); hold on

plot(data.distance,data.force_mx(:,2),'b')
plot(data.distance,data.force_vw(:,2),'r')
legend('maxwell stress tensor','virtual work')
title('Force per unit length y vs Distance x')
xlabel('Distance x (m)')
ylabel('Force per unit length (N/m)')


figure(6); hold on

plot(data.distance,data.torque_mx,'b')
plot(data.distance,data.torque_vw,'r')
legend('maxwell stress tensor','virtual work')
title('Torque per unit length z vs distance x')
xlabel('Distance x (m)')
ylabel('Torque per unit length (N*m/m)')

