%% Phase Difference with Respect to Reference Signal
%
% At this point I can find the direction of an incoming plane wave after I
% determine the relative phase shifts of the recordings on each antenna.
% The method used for finding phase shifts was applied in an ideal
% circumstance and is not applicable to the real problem at hand. In this
% code, I attempt to find the phase shift of a reference signal and another
% signal of the same frequency. The difference between the two phases
% should give me a relative phase shift which I can then use with my
% original direction algorithm.

clear all;
close all;
clc;

run('Sample_input_signal');
load('Sample_Antenna_Input.mat');


Eref = zeros(1,length(t));
for m=1:length(t)
        %Sine wave. Makes phase be zero
        Eref(m) = exp(1i*(omega*t(m)-pi/2));
end

figure;
plot(t, real(Eref));
axis([0 t(end) -E_o E_o]);

phsPts = floor(length(t)/numPeriods);
phsPtMax = length(t)-phsPts;

phase = zeros(1,phsPtMax);
bkr = zeros(4,phsPtMax);

for i=1:phsPtMax
    phase(i) = PhaseShift(real(Eref(i:(i+phsPts-1))),t(i:(i+phsPts-1))...
        ,omega);
    for j=1:4
        bkr(j,i)= PhaseShift(real(E(j,i:(i+phsPts-1))),t(i:(i+phsPts-1))...
            ,omega);
    end
end

figure;
hold on;
plot(phase);
hold on
plot(bkr(1,:));
hold on
plot(bkr(2,:));
hold on
plot(bkr(3,:));
hold on
plot(bkr(4,:));

