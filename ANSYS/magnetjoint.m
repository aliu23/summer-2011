%% 2D analysis of rotating magnet joint
%
% This file reads in ANSYS data and plots forces and torques vs angle for 2
% permanent magnets (mesh size 0.003).

%% Setup

clc

data=readansys('torqvsangle.ansys');

%%

figure(1); hold on

plot(data.angles,data.force_mx(:,1),'b')
plot(data.angles,data.force_vw(:,1),'r')
legend('maxwell stress tensor','virtual work')
title('Force per unit length x vs angle')
xlabel('Angle in degrees')
ylabel('Force per unit length (N/m)')


figure(2); hold on

plot(data.angles,data.force_mx(:,2),'b')
plot(data.angles,data.force_vw(:,2),'r')
legend('maxwell stress tensor','virtual work')
title('Force per unit length y vs angle')
xlabel('Angle in degrees')
ylabel('Force per unit length (N/m)')


figure(3); hold on

plot(data.angles,data.torque_mx,'b')
plot(data.angles,data.torque_vw,'r')
legend('maxwell stress tensor','virtual work')
title('Torque per unit length z vs angle')
xlabel('Angle in degrees')
ylabel('Torque per unit length (N*m/m)')


