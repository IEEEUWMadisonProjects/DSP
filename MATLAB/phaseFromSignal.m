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

%%

run('Sample_input_signal');
load('Sample_Antenna_Input.mat');
close all;

Eref = zeros(1,numTpoints);
for m=1:numTpoints
        %Sine wave. Makes phase be zero
        Eref(m) = exp(1i*(omega*t(m)-pi/2));
end

% figure;
% plot(t, real(Eref));
% axis([0 t(end) -E_o E_o]);

phsPts = floor(numTpoints/numPeriods);
phsPtMax = numTpoints-phsPts;

phase = zeros(1,phsPtMax);
bkr = zeros(4,phsPtMax);

for i=1:phsPtMax
    phase(i) = PhaseShift(real(Eref(i:(i+phsPts-1))),t(1:(phsPts))...
        ,omega);
    for j=1:4
        bkr(j,i)= PhaseShift(real(E(j,i:(i+phsPts-1))),t(1:(phsPts))...
            ,omega);
    end
end

phase = unwrap(phase);
for i=1:4
    bkr(i,:) = unwrap(bkr(i,:));
end

% figure;
% hold on;
% plot(phase);
% hold on
% plot(bkr(1,:));
% hold on
% plot(bkr(2,:));
% hold on
% plot(bkr(3,:));
% hold on
% plot(bkr(4,:));

figure;
plot(bkr(1,:)-phase);


%% Differece between phase section
% So now I have a running phase calculator and I want to apply it to a
% signal that is more realistic

% Ecomb = zeros(1,length(E)*4);
Ecomb = [E(1,:),E(2,:),E(3,:),E(4,:)];

figure;
plot(real(Ecomb));
axis([0 length(Ecomb) -4*E_o 4*E_o]);

t2 = 0:(t(2)-t(1)):numPeriods*T*4;

Eref2 = zeros(1,length(t2));
for m=1:length(t2)
        %Sine wave. Makes phase be zero
        Eref2(m) = exp(1i*(omega*t2(m)-pi/2));
end

hold on;
plot(t2, real(Eref2));
%axis([0 t2(end) -E_o E_o]);

phsPts = floor(length(t2)/(numPeriods*4));
phsPtMax = length(t2)-phsPts;

phase = zeros(1,phsPtMax);
bkr2 = zeros(1,phsPtMax);

for i=1:phsPtMax
    phase(i) = PhaseShift(real(Eref2(i:(i+phsPts-1))),t2(1:(phsPts))...
        ,omega);
    bkr2(i)= PhaseShift(real(Ecomb(i:(i+phsPts-1))),t2(1:(phsPts))...
           ,omega);
end

phase = unwrap(phase);
bkr2 = unwrap(bkr2);

figure;
plot(phase);
hold on;
plot(bkr2,'g');

figure;
plot(bkr2-phase);