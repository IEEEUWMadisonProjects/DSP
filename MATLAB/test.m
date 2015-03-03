%% Test file

%% Testing the PhaseShift Function

%Ensures we have the correct dataset to call the PhaseShift function
if(exist('Sample_Antenna_Input.mat','file')==0)
    run('Sample_input_signal');
end

%close all
%clear
%load relevant data
load('Sample_Antenna_Input.mat');

% find Phase shift between two sample signals
% First signal is reference signal
kr = zeros(1,4);
% kr(1) = 0;
% kr(2) = PhaseShift(real(E(1,:)), real(E(2,:)), t, omega);
% kr(3) = PhaseShift(real(E(1,:)), real(E(3,:)), t, omega);
% kr(4) = PhaseShift(real(E(1,:)), real(E(4,:)), t, omega);

kr(1) = PhaseShift(real(E(1,:)), t, omega);
kr(2) = PhaseShift(real(E(2,:)), t, omega);
kr(3) = PhaseShift(real(E(3,:)), t, omega);
kr(4) = PhaseShift(real(E(4,:)), t, omega);
%The kr values we have now are actually beta*k*r, so we need to divide by
%beta
kr = kr/beta;

knew = lsqlin(r_all,kr);
knew = knew/norm(knew);

figure;
quiver(0,0,knew(1),knew(2));
hold on;
scatter(r_all(:,1), r_all(:,2));
%axis([-.05 1.05 -.05 1.05]);
title('Guessed Direction');
xlabel('x');
ylabel('y');