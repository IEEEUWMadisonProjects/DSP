%% Test file 2
% Author: Alex Gabourie

%% Testing Finding Direction

close all


%Ensures we have the correct dataset to call the PhaseShift function
%if(exist('Sample_Antenna_Input.mat','file')==0)
    run('Sample_input_signal2');
%end

%close all
%clear
%load relevant data
load('Sample_Antenna_Input2.mat');

% find Phase shift between two sample signals
% First signal is reference signal
bkr = zeros(1,4);
bkr(1) = PhaseShift(real(E(1,:)), t, omega);
bkr(2) = PhaseShift(real(E(2,:)), t, omega);
bkr(3) = PhaseShift(real(E(3,:)), t, omega);
bkr(4) = PhaseShift(real(E(4,:)), t, omega);

%put the phases in the correct order and set the phase for first antenna to
%be 0 for the system of equations to be solved.
bkr = OrderPhase(bkr);

%The bkr values we have now are actually beta*k*r, so we need to divide by
%beta
kr = -bkr/beta;
r_n = r_all'*r_all;
knew = r_n\(r_all'*kr');
knew = knew/norm(knew);

figure;
quiver(0,0,knew(1),knew(2));
hold on;
scatter(r_all(:,1), r_all(:,2));
title('Ideal Signal Guessed Direction');
xlabel('x');
ylabel('y');

%% Second Stage Testing
%
% In this section, I will append the correct signals together in the way
% that our switching circuit will when running the real code. This means
% taking 2ms clips from each antenna in the order of 1,2,3,4. On the plot
% of the antenna locations this is starting at the origin and working
% counter clockwise. My plan is to use the 2ms to determine how much of a
% phase shift I would expect for the other antennas to have when I switch
% to them.

tVecIdx = zeros(1,4);

%not best code, just want something to work
%antenna 1 to 2 transition
for i=1:length(t)
    if t(i) > .002
        tVecIdx(1) = i;
        break;
    end
end
%antenna 2 to 3 transition
for i=tVecIdx(1):length(t)
    if t(i) > .004
        tVecIdx(2) = i;
        break;
    end
end   
%antenna 3 to 4 transition
for i=tVecIdx(2):length(t)
    if t(i) > .006
        tVecIdx(3) = i;
        break;
    end
end   
%antenna 4 to end transition
for i=tVecIdx(3):length(t)
    if t(i) > .008
        tVecIdx(4) = i;
        break;
    end
end     
    
    
    
    
    
    