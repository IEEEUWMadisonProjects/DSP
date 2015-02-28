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

%% Initialize Workspace
clear;
close all;
clc;

%% Input Parameters
%Frequency
f = 150e6; %[hz]

%E-field magnitude
E_o = 1; %[V/m] Not used as of now

%Vacuum permittivity 
epsilon_o = 8.8541878176e-12; %[F/m]
epsilon_rel = 1;

%permeability
mu_o = (4*pi)*10^(-7); %[H/m]
mu_rel = 1;
%recall that [F/m]*[H/m] = s^2/m^2 = 1/c

%The four antennas will be placed in a square grid like:
%  4    3
%
%  1    2
%
% The separation constant (like lattice constant) is 
a = 1; %[m]

%From this the positions of the antennas can be determined
r1 = [0, 0];
r2 = [a, 0];
r3 = [a, a];
r4 = [0, a];

%For plotting
r_all = [r1;r2;r3;r4];

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

%want 5 periods of information collected
T = 1/f;
t = linspace(0,5*T, 200); %time vector based on frequency

%Electric field matrix initialized to num time steps & antennas
E = zeros(length(bkr),length(t));

for m=1:length(t)
    for n=1:length(bkr)
        E(n,m) = exp(1i*(omega*t(m)-bkr(n)));
    end
end


numPlots = length(bkr);
gridSize = numPlots;
if sqrt(numPlots)==floor(sqrt(numPlots))
    gridSize = sqrt(numPlots);
else
    gridSize = floor(sqrt(numPlots))+1;
end
    
    
figure;
quiver(0,0,k(1),k(2));
hold on;
scatter(r_all(:,1), r_all(:,2));
axis([-.05 1.05 -.05 1.05]);
title('Initial Setup');
xlabel('x');
ylabel('y');

figure;
for m=1:length(bkr)
    subplot(gridSize,gridSize,m);
    plot(t,real(E(m,:)));
    title(['Antenna ' num2str(m)]);
    xlabel('Time [S]');
    ylabel('Amplitude [V/m]');
    axis([0 t(end) -E_o E_o]);
end

% I might want to end up saving the resulting information to a data file so
% that other parts of the DSP program can try to use the idealized output
% form this file and test their algorithms.

%% Plotting the Entire Plane Wave (t = 0)
% This section will create a two dimensional grid and then use the
% information selected above to create a three dimensional plot of the EM
% wave


%x and y coordinates
x = linspace(0,2*a,100);
y = linspace(0,2*a,100);

%new efield mesh to be calculated
E2D = zeros(length(y), length(x));

F(length(t)) = struct('cdata',[],'colormap',[]);
h = figure;

for p=1:length(t)
    for m=1:length(x)
        for n=1:length(y)
           xy = [x(m),y(n)];
           E2D(n,m) =  exp(1i*omega*t(p)-1i*beta*dot(k,xy));
        end 
    end
    
    surf(x,y,real(E2D))
    drawnow
    set(h,'Renderer','zbuffer') %workaround for bug. Not sure why it works
    F(p) = getframe(h);
end

writerObj = VideoWriter('PlaneWave.mp4', 'MPEG-4');
open(writerObj);
writeVideo(writerObj, F);
close(writerObj);