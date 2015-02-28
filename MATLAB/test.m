%% Test file

%% Testing the PhaseShift Function

%Ensures we have the correct dataset to call the PhaseShift function
if(exist('Sample_Antenna_Input.mat','file')==0)
    run('Sample_input_signal');
end

close all
clear
%load relevant data
load('Sample_Antenna_Input.mat');