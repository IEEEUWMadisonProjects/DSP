function [ ampPhaseOrdOut ] = OrderPhase2( ampPhaseOrd )
    for i=1:(length(ampPhaseOrd)-1)
        if ampPhaseOrd(3,i) > ampPhaseOrd(3,4)
            ampPhaseOrd(3,i) = ampPhaseOrd(3,i)-2*pi;
        end
    end
    ampPhaseOrdOut = ampPhaseOrd;
end

