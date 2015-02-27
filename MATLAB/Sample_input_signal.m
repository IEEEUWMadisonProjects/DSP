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

%% Input Parameters
%Frequency
f = 150e6; %[hz]

%E-field magnitude
E_o = 1; %[V/m]

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
a = .5; %[m]

%From this the positions of the antennas can be determined
r1 = [0, 0];
r2 = [a, 0];
r3 = [a, a];
r4 = [0, a];

%Propagation direction (theta = 0 => dirction is from antenna 1 to 2)
theta = 60; %[degrees]


%% Calculations

%angular frequency based off frequency
omega = 2*pi*f;

%permitivity & permeability (same units)
eps = epsilon_o*epsilon_rel;
mu = mu_o*mu_rel;

%propagation constant
beta = omega*sqrt(mu*eps); %[m^-1]

% The signal that we will be looking at is a 2D E-field of an EM wave
% propagating with no phase shift between H & E.


