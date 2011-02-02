%% 2D analysis of rotating magnet joint
%
% This file reads in ANSYS data and plots forces and torques vs angle of 2
% permanent magnets (mesh size 0.003).

%% Setup

clc

data=readansys('torqvsangle.ansys');

%%

figure(1); hold on

plot(data.angles,data.force_mx(:,1),'b')
plot(data.angles,data.force_vw(:,1),'r')
legend('maxwell stress tensor','virtual work')
title('Force x vs angle')
xlabel('Angle in degrees')
ylabel('Force in Newtons')


figure(2); hold on

plot(data.angles,data.force_mx(:,2),'b')
plot(data.angles,data.force_vw(:,2),'r')
legend('maxwell stress tensor','virtual work')
title('Force Y vs angle')
xlabel('Angle in degrees')
ylabel('Force in Newtons')


figure(3)

plot(data.angles,data.torque_mx,'b',data.angles,data.torque_vw,'r')
legend('maxwell stress tensor','virtual work')
title('Torque Z vs angle')
xlabel('Angle in degrees')
ylabel('Torque in Newton metres')


