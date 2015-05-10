%% Test 3
% playing with fabricated signals and moving antennas.

close all
clear

omega = 2*pi*150e6; %[rads/s]

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
a = .4699*sqrt(2); %[m] - about what our antenna spacing is
%From this the positions of the antennas can be determined
% shiftx = 0;%7*0.707596301024091;
% shifty = 0;%7*0.707596301024091;
% shift = [shiftx,shifty];
% r1 = [0, 0]+shift;
% r2 = [a, 0]+shift;
% r3 = [a, a]+shift;
% r4 = [0, a]+shift;
r1 = [-a/2, -a/2];
r2 = [a/2, -a/2];
r3 = [a/2, a/2];
r4 = [-a/2, a/2];

%For plotting
r_all = [r1;r2;r3;r4];
 
%% E

%Propagation direction (theta = 0 => dirction is from antenna 1 to 2 and 
%   increasing theta moves counter-clockwise)
theta = 2*pi*270/360; %[radians]
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
T = 2*pi/omega;
numTpoints = 48e3*.02; % samples/sec*seconds
t = linspace(0,.02,numTpoints); %time vector based on frequency [s]

%Electric field matrix initialized to num time steps & antennas
superDuperPhaseShift = 3;
E = zeros(length(bkr),length(t));
alpha = 1;
for m=1:length(t)
    for n=1:length(bkr)
        E(n,m) = imag(exp(1i*(omega*t(m)-bkr(n)+superDuperPhaseShift))*exp(-bkr(n)/beta*alpha));
    end
end

center = r3-(r3-r1)/2;
figure;
quiver(center(1),center(2),-k(1),-k(2));
hold on;
scatter(r_all(:,1), r_all(:,2));
title('Initial Setup');
xlabel('x');
ylabel('y');
axis([center(1)-1.5 center(1)+1.5 center(2)-1.5 center(2)+1.5]);

%% Ein    
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

%% Eref

 Eref = zeros(1,length(Ein));
    for m=1:length(Ein)
            %Sine wave. Makes phase be zero
            Eref(m) = real(exp(1i*(omega*t(m)-pi/2)));
    end

    
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
    
    % For proof of concept. If number of points is large, this should
    % be close enough
    phsAveSize = floor(phsPtMax/8); 
    %The number of points put in to the average that gives us our phase
    % for each antenna

    startIndex = zeros(1,4);
    for i=1:4
        startIndex(i) = floor(phsPtMax/4)*i-floor(phsAveSize*1.5);
    end
    % Since the antennas are switched from 4,3,2,1 we need to reverse the
    % order of our start indexes  and then use those in the bkr stuff.
    bkr = [mean(phaseDiff(startIndex(1):(startIndex(1)+phsAveSize))),...
        mean(phaseDiff(startIndex(2):(startIndex(2)+phsAveSize))),...
        mean(phaseDiff(startIndex(3):(startIndex(3)+phsAveSize))),...
        mean(phaseDiff(startIndex(4):(startIndex(4)+phsAveSize)))];
    bkr/pi
    OrderPhase(bkr/pi)
    if ~(sum(bkr)<10^(-4) && sum(bkr)>-10^(-4))
        bkr = bkr - sum(bkr)/4;
    end
    
    %bkr = bkr - bkr(1)

    amp = [mean(expAmp(startIndex(1):(startIndex(1)+phsAveSize))),...
        mean(expAmp(startIndex(2):(startIndex(2)+phsAveSize))),...
        mean(expAmp(startIndex(3):(startIndex(3)+phsAveSize))),...
        mean(expAmp(startIndex(4):(startIndex(4)+phsAveSize)))];

    %possible refactor needed
    ampPhaseOrd = sortrows([amp;linspace(1,length(amp),length(amp));bkr]',1)';
    ampPhaseOrd(3,:)/pi
    ampPhaseOrd(3,:) = OrderPhase(ampPhaseOrd(3,:));
    ampPhaseOrd(3,:)/pi
    
        % Take the first two indicies from ampPhaseOrd to determine a direction
    % I know if directly on:
    % 1: 45
    % 2: 135
    % 3: 225
    % 4: 315
    % Rotations are based on angle from vector drawn from antenna 1 to 2
    % Recognizing which antenna has the largest amplitude narrows it to a 45
    % degree angle
    center = 45+90*(ampPhaseOrd(2,4)-1);
    % with the center angle, you either add up to 45 degrees or subtract up to
    % 45 degrees
    antNum = [1,2,3,4,1]; %Didn't feel like making circular buffer
    sign = 1;

    if (ampPhaseOrd(2,4)==1 && ampPhaseOrd(2,3)==4) || ...
       (ampPhaseOrd(2,4)==2 && ampPhaseOrd(2,3)==1) || ...
       (ampPhaseOrd(2,4)==3 && ampPhaseOrd(2,3)==2) || ...
       (ampPhaseOrd(2,4)==4 && ampPhaseOrd(2,3)==3)

        sign=-1;
    end
    topSig = ampPhaseOrd(1,4)-ampPhaseOrd(1,2);
    secondSig = ampPhaseOrd(1,3)-ampPhaseOrd(1,2);
    inAngle = (center+90*sign*(secondSig/(topSig+secondSig)));
    inAngle = inAngle*pi/180;
    knew3 = [1,0];
    rot = [cos(inAngle),-sin(inAngle);sin(inAngle),cos(inAngle)];
    knew3 = -(rot*knew3')';
    
    
    %%%%% 3 Ant Method %%%%%
    % Use the three strongest antennas, in strongest to weakest order
    rcurr = [r_all(ampPhaseOrd(2,4),:);
             r_all(ampPhaseOrd(2,3),:);
             r_all(ampPhaseOrd(2,2),:)];

    %using strongest to weakest signal strength, order the phases
    bkrTop3 = [ampPhaseOrd(3,4),ampPhaseOrd(3,3),ampPhaseOrd(3,2)];
    krTop3 = bkrTop3/beta;
    r_nTop3 = rcurr'*rcurr;
    knew4 = r_nTop3\(rcurr'*krTop3');
    knew4 = knew4/norm(knew4);
    
    %%%%% 2 Ant Method %%%%%
    % Use the two strongest antennas, in strongest to weakest order
    rcurr = [r_all(ampPhaseOrd(2,4),:);
             r_all(ampPhaseOrd(2,3),:)];

    %using strongest to weakest signal strength, order the phases
    bkrTop2 = [ampPhaseOrd(3,4),ampPhaseOrd(3,3)];
    
    krTop2 = bkrTop2/beta;
    knew5 = rcurr\krTop2';
    knew5 = knew5/norm(knew5);
    
    
        %%%%% PHASE PLOT %%%%%
    figure1 = figure;
    % plot(t(1:length(phaseDiff))*1000, phaseDiff);
    plot(phaseDiff);
    title('Phase Difference vs. Time','FontSize',14);    
    ylabel('Phase Difference [Rads]','FontSize',12);
    xlabel('Time [ms]','FontSize',12);
    padding = .15*max(abs(min(phaseDiff)),abs(max(phaseDiff)));
    axis([0 phsPtMax min(phaseDiff)-padding max(phaseDiff)+padding ]);
    set(figure1, 'Position', [1000, 200, 800, 600]);
    
        %%%%% AMPLITUDE PLOT %%%%%
    figure2 = figure;
    plot(expAmp);
    title('Amplitude vs. Time','FontSize',14);    
    ylabel('Amplitude [Arb]','FontSize',12);
    xlabel('Time [ms]','FontSize',12);
    axis([0 phsPtMax 0 1.15*max(expAmp)]);
    set(figure2, 'Position', [100, 200, 800, 600]);
    
        %%%%% FINAL PLOT %%%%%
    %Get line segments for axis to relate to Tom's figures
    line13x = [r1(1) r3(1)];
    line13y = [r1(2) r3(2)];
    line24x = [r2(1) r4(1)];
    line24y = [r2(2) r4(2)];
    
    figure3 = figure;
    center = r3-(r3-r1)/2;
    % Real k
    k= -k;
    quiver(center(1),center(2), k(1), k(2));
    hold on
%     % 4 Antenna Phase Method
%     quiver(center(1),center(2),knew(1),knew(2));
%     hold on; 
%     % Attenuation Method
%     quiver(a/2,a/2,knew2(1),knew2(2));
%     hold on;
    % Amplitude Method - 3 Antennas
    quiver(center(1),center(2),knew3(1),knew3(2));
    hold on;
%     % 3 Antenna Phase Method (top 3 amplitudes)
    quiver(center(1),center(2),knew4(1),knew4(2));
    hold on;
    % 2 Antenna Phase Method (top 2 amplitudes)
    quiver(center(1),center(2),knew5(1),knew5(2));
    hold on;
    
    scatter(r_all(:,1), r_all(:,2));
%     title([fileName ': Guessed Direction '],'FontSize',14);
    xlabel('x [m]','FontSize',14);
    ylabel('y [m]','FontSize',14);
%     axis([-1 2 -1 2]);
    axis([center(1)-1.5 center(1)+1.5 center(2)-1.5 center(2)+1.5]);
    legend('on');
%     legend('4 Ant','Amplitude', 'Phase(3)', 'Phase(2)');
    legend('Real','Amplitude', 'Phase(3)', 'Phase(2)');
%     legend('Amplitude',  'Phase(2)');
    %Draw coordinate lines
    plot(line13x,line13y);
    hold on
    plot(line24x,line24y);
    hold on
    
    createAntennaTextBox(figure3);
    set(figure3, 'Position', [200, 200, 800, 800]);