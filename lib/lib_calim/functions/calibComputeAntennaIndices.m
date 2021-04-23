function [Pos,Posc,A1S2,A2S1] = calibComputeAntennaIndices(SuppD2,SuppD1,Antennas,flag,param)

nTimeSlots     = param.nTimeSlots;
InactiveTimeSlot = param.InactiveTimeSlot;
AntsIndices     = param.AntennasIndex ;
nAntennas     = param.nAntennas;

Pos  = cell(nAntennas,nTimeSlots);
Posc = cell(nAntennas,nTimeSlots);
A1S2 = cell(nAntennas,nTimeSlots);
A2S1 = cell(nAntennas,nTimeSlots);

% gather indices of the antennas
for iAnt = 1: nAntennas
    for itimeSlot = 1:nTimeSlots
        flagTemp =  flag{itimeSlot};
        if ~ismember(itimeSlot,InactiveTimeSlot)
            [~,A1S2{iAnt,itimeSlot}]=ismember(Antennas.ant1{itimeSlot}(SuppD2{iAnt,itimeSlot}),AntsIndices);
            A1S2{iAnt,itimeSlot} = (A1S2{iAnt,itimeSlot});
            if ~isempty(A1S2{iAnt,itimeSlot})
                Pos{iAnt,itimeSlot}=(flagTemp(SuppD2{iAnt,itimeSlot}));
            end
            [~,A2S1{iAnt,itimeSlot}]=ismember(Antennas.ant2{itimeSlot}(SuppD1{iAnt,itimeSlot}),AntsIndices);
            A2S1{iAnt,itimeSlot} = (A2S1{iAnt,itimeSlot});
            if ~isempty(A2S1{iAnt,itimeSlot})
                Posc{iAnt,itimeSlot}=(flagTemp(SuppD1{iAnt,itimeSlot}));
            end
        end            
      
    end
end
end
