%% Simple Model of Antenna Array
% The purpose of this script is to create a model of the signals that we 
% might receive from the radio location collars. We assume that, since the
% antennas are far enough away from the source, we see the incoming waves
% as plane waves. We also will only work in two dimensions, which is as if
% the waves and the antennas are on the same plane.
%
% Hopefully from this we can build a fake set of signals that allow us to
% understand how to interpret phase changes to find the direction of the
% wave propagation.

%Author: Alex Gabourie

%% EDIT
%
% This version is similar to the original file but some physical parameters
% will be changed as well as vector sizes. My goal for this is to better
% mimic the 20 ms recording times that we would expect for the actual
% incoming signal

%% Initialize Workspace
clear;
close all;
clc;

runVideo = 0;

%% Input Parameters
%Frequency
f = 2551; %[hz] - measured from FFT of signal from gonzo file

%E-field magnitude
E_o = 1; %[V/m] Not used as of now

%Vacuum permittivity 
epsilon_o = 8.8541878176e-12; %[F/m]
epsilon_rel = 1;

%permeability
mu_o = (4*pi)*10^(-7); %[H/m]
mu_rel = 1;
%recall that [F/m]*[H/m] = s^2/m^2 = 1/c

%Propagation direction (theta = 0 => dirction is from antenna 1 to 2 and 
%   increasing theta moves counter-clockwise)
theta = 2*pi*60/360; %[radians]

%% Calculations

%angular frequency based off frequency
omega = 2*pi*f;

%permitivity & permeability (same units)
eps = epsilon_o*epsilon_rel;
mu = mu_o*mu_rel;

%propagation constant
beta = omega*sqrt(mu*eps); %[m^-1]

%wavelength for intuition purposes
lambda = (2*pi)/beta;

%The four antennas will be placed in a square grid like:
%  4    3
%
%  1    2
%
% The separation constant (like lattice constant) is 
a = 1; %[m] - about what our antenna spacing is
%This constant is calculated to place the antennas an appropriate distance
%away such that we can determine which antenna was hit first

%From this the positions of the antennas can be determined
r1 = [0, 0];
r2 = [a, 0];
r3 = [a, a];
r4 = [0, a];

%For plotting
r_all = [r1;r2;r3;r4];

%Propagation direction vector
k = [1,0];
rot = [cos(theta),-sin(theta);sin(theta),cos(theta)];
k = (rot*k')';

%prop constant dotted with r as in time-harmonic exponential. We are only
% looking at a few positions, so these can be pre-calculated
bkr = beta*[dot(k,r1);
            dot(k,r2);
            dot(k,r3);
            dot(k,r4)];

%want 20ms of info collected
T = 1/f;
numTpoints = 48e3*.02; % samples/sec*seconds
t = linspace(0,.02,numTpoints); %time vector based on frequency [s]

%Electric field matrix initialized to num time steps & antennas
E = zeros(length(bkr),length(t));

for m=1:length(t)
    for n=1:length(bkr)
        E(n,m) = exp(1i*(omega*t(m)-bkr(n)));
    end
end


% numPlots = length(bkr);
% gridSize = numPlots;
% if sqrt(numPlots)==floor(sqrt(numPlots))
%     gridSize = sqrt(numPlots);
% else
%     gridSize = floor(sqrt(numPlots))+1;
% end
    
    
figure;
quiver(0,0,k(1),k(2));
hold on;
scatter(r_all(:,1), r_all(:,2));
title('Initial Setup');
xlabel('x');
ylabel('y');

% Plot works, no need for it anymore

% figure;
% for m=1:length(bkr)
%     subplot(gridSize,gridSize,m);
%     plot(t,real(E(m,:)));
%     title(['Antenna ' num2str(m)]);
%     xlabel('Time [S]');
%     ylabel('Amplitude [V/m]');
%     axis([0 t(end) -E_o E_o]);
% end

% Save real part of E-field to use in another script. This way we can
% interpret the information independent of the set of data in this script.
save('Sample_Antenna_Input2.mat', 'E', 'a', 'omega', 'r_all', 'k', 't', ...
        'beta', 'lambda');

