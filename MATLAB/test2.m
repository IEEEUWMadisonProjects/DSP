%% Test file 2
% Author: Alex Gabourie



close all
clear;
clc;

% input=0 -> sample input signal, input!=0 -> .wav
input = 1;

%% Testing Finding Direction
if(input == 0)
    %Ensures we have the correct dataset to call the sineFit function
    %if(exist('Sample_Antenna_Input.mat','file')==0)
        run('Sample_input_signal2');
    %end

    %close all
    %clear
    %load relevant data
    load('Sample_Antenna_Input2.mat');
    input = 0;

    % find Phase shift between two sample signals
    % First signal is reference signal
    bkr = zeros(1,4);
    bkr(1) = sineFit(real(E(1,:)), t, omega);
    bkr(2) = sineFit(real(E(2,:)), t, omega);
    bkr(3) = sineFit(real(E(3,:)), t, omega);
    bkr(4) = sineFit(real(E(4,:)), t, omega);

    %put the phases in the correct order and set the phase for first 
    %antenna to be 0 for the system of equations to be solved.
    bkr = OrderPhase(bkr);

    %The bkr values we have now are actually beta*k*r, so we need to divide
    %by beta
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

else
%     Ein = 25*audioread('fourAntennaSnip.wav');
%     Ein = audioread('1_36_1_Jason_Walk_single.wav');%2753
%    Ein = audioread('7_09_4_Jason_Walk_Single.wav'); %2683 Hz
%     Ein = audioread('10_12_3_Jason_Walk_Single.wav'); %2683

%%%%%%% WedApr22 %%%%%%%%%%%%%%%%
% Ein = audioread('Single_0deg.wav');omega = 2*pi*2612; %[rads/s]
% Ein = audioread('Single_30deg.wav');omega = 2*pi*2612; %[rads/s]
% Ein = audioread('Single_60deg.wav');omega = 2*pi*2663; %[rads/s]
% Ein = audioread('Single_90deg.wav');omega = 2*pi*2612; %[rads/s]
% Ein = audioread('Single_120deg.wav');omega = 2*pi*2700; %[rads/s]
% Ein = audioread('Single_150deg.wav');omega = 2*pi*2700; %[rads/s]
% Ein = audioread('Single_180deg.wav');omega = 2*pi*2715; %[rads/s]
Ein = audioread('Single_210deg.wav');omega = 2*pi*2715; %[rads/s]


%     omega = 2*pi*2612; %[rads/s]
    %want 20ms of info collected
    T = 2*pi/omega;
    numTpoints = 48e3*.02; % samples/sec*seconds
    t = linspace(0,.02,numTpoints); %time vector based on frequency [s]
    
    %Vacuum permittivity 
    epsilon_o = 8.8541878176e-12; %[F/m]
    epsilon_rel = 1;
    %permeability
    mu_o = (4*pi)*10^(-7); %[H/m]
    mu_rel = 1;
    
    %permitivity & permeability (same units)
    eps = epsilon_o*epsilon_rel;
    mu = mu_o*mu_rel;

    %propagation constant
    beta = omega*sqrt(mu*eps); %[m^-1]
    
    % The separation constant (like lattice constant) is 
    a = 1; %[m] - about what our antenna spacing is
    %From this the positions of the antennas can be determined
    r1 = [0, 0];
    r2 = [a, 0];
    r3 = [a, a];
    r4 = [0, a];

    %For plotting
    r_all = [r1;r2;r3;r4];

end
    
    
%% Second Stage Testing
%
% In this section, I will append the correct signals together in the way
% that our switching circuit will when running the real code. This means
% taking 2ms clips from each antenna in the order of 1,2,3,4. On the plot
% of the antenna locations this is starting at the origin and working
% counter clockwise. My plan is to use the 2ms to determine how much of a
% phase shift I would expect for the other antennas to have when I switch
% to them.

if(input == 0)
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

    % compose the antenna signal based on the time the antenna readings come in

    Ein = [real(E(1,1:tVecIdx(1)-1)), real(E(2,tVecIdx(1):tVecIdx(2)-1)),...
        real(E(3,tVecIdx(2):tVecIdx(3)-1)),...
        real(E(4,tVecIdx(3):tVecIdx(4)-1))];
end
    
Eref = zeros(1,length(Ein));
for m=1:length(Ein)
        %Sine wave. Makes phase be zero
        Eref(m) = real(exp(1i*(omega*t(m)-pi/2)));
end

%Plot the signals coming in from the antenna vs the reference signal
% figure;
% plot(Ein,'b');
% hold on;
% plot(Eref,'r');
% axis([0 length(Eref) -4 4]);
    

phsPts = floor(length(Eref)/(5*4));
phsPtMax = length(Eref)-phsPts;

%phase from reference sine wave
realPhase = zeros(1,phsPtMax);
% experimental phase
expPhase = zeros(1,phsPtMax);
expAmp = zeros(1,phsPtMax);%to get amplitude plot much like the phase


for i=1:phsPtMax
    realPhase(i) = sineFit(real(Eref(i:(i+phsPts-1))),t(1:(phsPts))...
        ,omega);
    [expPhase(i),expAmp(i)] = sineFit(real(Ein(i:(i+phsPts-1))),...
        t(1:(phsPts)),omega);
end

%unwrap phases to get the phase difference plot
realPhase = unwrap(realPhase);
expPhase = unwrap(expPhase);

% Phase differences
phaseDiff = expPhase-realPhase;


figure1 = figure;
% plot(t(1:length(phaseDiff))*1000, phaseDiff);
plot(phaseDiff);
title('Phase Difference vs. Time','FontSize',14);    
ylabel('Phase Difference [Rads]','FontSize',12);
xlabel('Time [ms]','FontSize',12);
% axis([0 phsPtMax -3*pi 3*pi]);
% axis([0 phsPtMax 2 4]);
% axis([0 t(length(phaseDiff))*1000 0 1.5*pi]);
padding = .15*max(abs(min(phaseDiff)),abs(max(phaseDiff)));

axis([0 phsPtMax min(phaseDiff)-padding max(phaseDiff)+padding ]);
% annotation(figure1,'textbox',...
%     [0.174214285714286 0.797619047619048 0.190071428571429 0.0666666666666721],...
%     'String',{'freq. = 2753 Hz'},...
%     'FitBoxToText','off',...
%     'EdgeColor',[0.941176474094391 0.941176474094391 0.941176474094391]);

figure2 = figure;
plot(expAmp);
title('Amplitude vs. Time','FontSize',14);    
ylabel('Amplitude [Arb]','FontSize',12);
xlabel('Time [ms]','FontSize',12);
axis([0 phsPtMax 0 1.15*max(expAmp)]);


%%%%%%%%%%%%%%%%%%% Direction commented out until we determine if the phase
%%%%%%%%%%%%%%%%%%% method even makes sense
% % At this point we should see the phase differences that correspond to a
% % particular direction and that resembles reality to a decent degree. Since
% % this is an idealized situation, our phase difference plots have some
% % flat regions from which we can sample from.
% 
% 
% % Obviously this needs to be done differently, but it is for proof of
% % concept
% % bkr = [phaseDiff(32), phaseDiff(141), phaseDiff(235), phaseDiff(335)];
% % bkr = [phaseDiff(26), phaseDiff(82), phaseDiff(160), phaseDiff(216)];
% bkr = [phaseDiff(26), phaseDiff(82), phaseDiff(144), phaseDiff(203)];
% 
% %put the phases in the correct order and set the phase for first antenna to
% %be 0 for the system of equations to be solved.
% bkr = OrderPhase(bkr);
% 
% %The bkr values we have now are actually beta*k*r, so we need to divide by
% %beta
% kr = -bkr/beta;
% r_n = r_all'*r_all;
% knew = r_n\(r_all'*kr');
% knew = knew/norm(knew);
% 
% figure;
% quiver(0,0,knew(1),knew(2));
% hold on;
% scatter(r_all(:,1), r_all(:,2));
% title('Guessed Direction','FontSize',14);
% xlabel('x [m]','FontSize',14);
% ylabel('y [m]','FontSize',14);
% axis([-.6 1.1 -1 1.1]);


%% Notes
% I need to discuss with Tom how the signals are being recorded because I
% still have a little confusion. Its the demodulation that trips me up
% because in the 20 ms timeframe I would expect to have larger phase
% differences between the signals of each antenna. Since the signal is
% demodulated to ~2551 Hz the phase differences are miniscule, but luckily,
% with double precision numbers, enough to capture in this simulation.
% However, this worries me for real testing where small inaccuracies may
% mess with the entire phase difference determinations
    