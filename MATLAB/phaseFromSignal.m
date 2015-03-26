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

run('Sample_input_signal');
load('Sample_Antenna_Input.mat');

