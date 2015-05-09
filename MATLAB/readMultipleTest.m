%% Read Multiple Files Test
%
% Author: Alex Gabourie
% Reading multiple files and test their direction
%% Initialize workspace
close all
clear;
clc;

%% User input parameters
%The folder in which all of the wav files are found
folderNum = 0;

if folderNum==0
    folderName = '.\May_1_4Ant';
    throwAwayPoints = 50;
elseif folderNum==1
    folderName = '.\apr22';
    throwAwayPoints = 0;
end

folderOut = folderName;

%number of antennas used in signal
numAnt = 4;

%Grabs the file list and places them into the filesList variable
filesList = dir([folderName '\*.wav']);

%wav files are expected to contain only the information of interest i.e.
%only the readings from the number of antennas during their cycle. It may
%be the case that there are a particular number of points before the first
%antenna reading. In this case list the number of points to throw out from
%the wav file in the beginning



%% Constant parameters
    
    %Vacuum permittivity 
    epsilon_o = 8.8541878176e-12; %[F/m]
    epsilon_rel = 1;
    %permeability
    mu_o = (4*pi)*10^(-7); %[H/m]
    mu_rel = 1;
    
    %permitivity & permeability (same units)
    eps = epsilon_o*epsilon_rel;
    mu = mu_o*mu_rel;
    
    % The separation constant (like lattice constant) is 
    a = 1; %[m] - about what our antenna spacing is
    %From this the positions of the antennas can be determined
    r1 = [0, 0];
    r2 = [a, 0];
    r3 = [a, a];
    r4 = [0, a];

    %For plotting
    r_all = [r1;r2;r3;r4];

%% File Loop
    
for fileNum=5:5%length(filesList)
    close all;
    fileName = filesList(fileNum).name;
    file = [folderName '\' fileName];
    Ein = audioread(file); 
    Ein = Ein(throwAwayPoints+1:end,1)';
    
    % The DTFT (FFT) yields a continuous spectrum of frequencies that we
    % must sample to get a finite, representable number of frequencies in
    % the computer. However, you can choose the number of points you want
    % in your FFT to get whatever resolution you would like in the
    % frequency domain.
    NFFT = 2^(nextpow2(length(Ein))+2);
    EinFFT = fft(Ein(1:floor(length(Ein)/4)),NFFT);
    [val, idx] = max(abs(EinFFT(1:floor(NFFT/2))));
    freq = idx*2*pi/(NFFT*2*pi)*48000;
    freqString = num2str(floor(freq))
%     figure
%     plot(abs(EinFFT(1:floor(length(EinFFT)/2))));
    omega = 2*pi*(freq);
    
    %want 20ms of info collected
    T = 2*pi/omega;
    numTpoints = 48e3*.02; % samples/sec*seconds
    t = linspace(0,.02,numTpoints); %time vector based on frequency [s]
    %propagation constant
    beta = omega*sqrt(mu*eps); %[m^-1]
    
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
    bkr = [mean(phaseDiff(startIndex(4):(startIndex(4)+phsAveSize))),...
        mean(phaseDiff(startIndex(3):(startIndex(3)+phsAveSize))),...
        mean(phaseDiff(startIndex(2):(startIndex(2)+phsAveSize))),...
        mean(phaseDiff(startIndex(1):(startIndex(1)+phsAveSize)))];


    amp = [mean(expAmp(startIndex(4):(startIndex(4)+phsAveSize))),...
        mean(expAmp(startIndex(3):(startIndex(3)+phsAveSize))),...
        mean(expAmp(startIndex(2):(startIndex(2)+phsAveSize))),...
        mean(expAmp(startIndex(1):(startIndex(1)+phsAveSize)))];

    %possible refactor needed
    ampPhaseOrd = sortrows([amp;linspace(1,length(amp),length(amp));bkr]',1)';

    
    %%%%% 4 Ant Phase Method %%%%
    % %The bkr values we have now are actually beta*k*r, so we need to divide by
    % %beta
    % kr = -bkr/beta;
    % r_n = r_all'*r_all;
    % knew = r_n\(r_all'*kr');
    % knew = knew/norm(knew);
    
    %%%%% ATTENUATION METHOD %%%%%
    %%%%% DIDN'T WORK %%%%%
    % %Find direction through using the amplitides
    % alpha = 1.9521; %an empirically found attenuation constant
    % %B*exp(-alpha*dot(k,r)) = C where C is amp
    % B = max(amp);
    % linAmp = log(amp/B)./(-alpha);
    % % linAmp = amp;
    % r_n = r_all'*r_all;
    % knew2 = r_n\(r_all'*linAmp');
    % knew2 = -knew2/norm(knew2);
    
    
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
    % Use the three strongest antennas, in strongest to weakest order
    rcurr = [r_all(ampPhaseOrd(2,4),:);
             r_all(ampPhaseOrd(2,3),:)];

    %using strongest to weakest signal strength, order the phases
    bkrTop2 = [ampPhaseOrd(3,4),ampPhaseOrd(3,3)];
    krTop2 = -bkrTop2/beta;
    r_nTop2 = rcurr'*rcurr;
    knew5 = r_nTop2\(rcurr'*krTop2');
    knew5 = knew5/norm(knew5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %%%%% REF VS EIN %%%%%
%     figure0 = figure;
%     plot(Ein,'b');
%     hold on;
%     plot(Eref,'r');
%     axis([0 length(Eref) -4 4]);
    
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
%     annotation(figure1,'textbox',...
%         [0.174214285714286 0.797619047619048 0.190071428571429 ...
%         0.0666666666666721],'String',{['f = ' freqString ' Hz']},...
%         'FitBoxToText','off',...
%         'EdgeColor',[0.941176474094391 0.941176474094391 ...
%         0.941176474094391]);
    
    %%%%% AMPLITUDE PLOT %%%%%
%     figure2 = figure;
%     plot(expAmp);
%     title('Amplitude vs. Time','FontSize',14);    
%     ylabel('Amplitude [Arb]','FontSize',12);
%     xlabel('Time [ms]','FontSize',12);
%     axis([0 phsPtMax 0 1.15*max(expAmp)]);
%     set(figure2, 'Position', [100, 200, 800, 600]);
    
    %%%%% FINAL PLOT %%%%%
    %Get line segments for axis to relate to Tom's figures
    line13x = [r1(1) r3(1)];
    line13y = [r1(2) r3(2)];
    line24x = [r2(1) r4(1)];
    line24y = [r2(2) r4(2)];
    
    figure3 = figure;
%     % 4 Antenna Phase Method
%     quiver(a/2,a/2,knew(1),knew(2));
%     hold on; 
%     % Attenuation Method
%     quiver(a/2,a/2,knew2(1),knew2(2));
%     hold on;
    % Amplitude Method - 3 Antennas
    quiver(a/2,a/2,knew3(1),knew3(2));
    hold on;
    % 3 Antenna Phase Method (top 3 amplitudes)
    quiver(a/2,a/2,knew4(1),knew4(2));
    hold on;
    % 2 Antenna Phase Method (top 2 amplitudes)
    quiver(a/2,a/2,knew5(1),knew5(2));
    hold on;
    
    scatter(r_all(:,1), r_all(:,2));
    title([fileName ': Guessed Direction '],'FontSize',14);
    xlabel('x [m]','FontSize',14);
    ylabel('y [m]','FontSize',14);
    axis([-1 2 -1 2]);
    legend('on');
    legend('Amplitude', 'Phase(3)', 'Phase(2)');
    %Draw coordinate lines
    plot(line13x,line13y);
    hold on
    plot(line24x,line24y);
    hold on
    
    createAntennaTextBox(figure3);
    set(figure3, 'Position', [200, 200, 800, 800]);
    fileTemp = strsplit(fileName, '.wav');
    if ~(exist([folderOut '_out'],'dir')==7)
        mkdir([folderOut '_out']);
    end
    fileTemp = regexprep(strjoin([folderOut '_out\' fileTemp(1) '_dir.png']),'\s','');
    saveas(figure3,fileTemp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end