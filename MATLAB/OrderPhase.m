function [ phaseRet ] = OrderPhase( phases )
%ORDERPHASE Fixes the ordering of phase information from antennas. 
%   Takes the phases from the antenna signals and determines if
%   they are in the correct order. If not, it fixes the ordering. We are
%   looking for differences between two phases that exceed pi/2 (quarter
%   wavelength)

%ordered matrix
ord = sortrows([phases;linspace(1,length(phases),length(phases))]',1)';

indx = -1;

%move through ordered array looking for big differences in phase
for i=1:length(phases)-1
    % if there is a difference in phase, make note and store which antennas
    if (ord(1,i+1) - ord(1,i)) >= pi
        indx = i;
        break;
    end
end

% force the signals detected to be first to have the highest phase
if indx > 0
    for i=indx:-1:1
        ord(1,i) = ord(1,i)+2*pi;
    end
end

%sort the data to its origninal order and return the new phases
ord = (sortrows(ord',2)');
if(max(ord(1,:))>pi)
    ord(1,:) = ord(1,:) - (max(ord(1,:))-pi);
end
phaseRet = ord(1,:);
