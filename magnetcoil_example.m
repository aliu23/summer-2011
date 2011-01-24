%% Force between a coil and a magnet
%
% First of all, a magnet inside a coil:

close all
clear all
clc

%% Test for correctness (I hope!)

f = magnetcoil(10e-3,9e-3,10e-3,1,10e-3,10.5e-3,20e-3,100,1);
assert( round(100000*f)==-62946 , 'Test case failed!')

%% Magnet parameters

rm = 0.009; % radius
lm = 0.010; % length
Br = 1.3;   % magnet strength for a rare-earth magnet


%% Coil parameters

rc = 0.01;  % inner radius
Rc = 0.015; % outer radius
lc = 0.02;  % length

N = 100; % turns of wire
I = 1;   % current (amps)

%% Calculations and plot

% axial displacement range between magnet/coil centres:
z = linspace(0.001,0.04);

Fz = magnetcoil(z,rm,lm,Br,rc,Rc,lc,N,I); % calculate forces

plot(z,Fz)


%% Force between a pancake coil and a magnet
%
% First of all, a magnet next to a pancake coil:

%% Magnet parameters

rm = 0.02; % radius
lm = 0.01; % length
Br = 1.3;  % magnet strength for a rare-earth magnet


%% Coil parameters

rc = 0.005; % inner radius
Rc = 0.015; % outer radius
lc = 0.01;  % length

N = 100; % turns of wire
I = 1;   % current (amps)

%% Calculations and plot

% axial displacement range between magnet/coil faces:
z = lm/2+lc/2+linspace(0.001,0.04);

Fz = magnetcoil(z,rm,lm,Br,rc,Rc,lc,N,I); % calculate forces

figure
plot(z,Fz)