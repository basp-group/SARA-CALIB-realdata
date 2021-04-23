function [data_spec,gVars,snap_gen] = util_get_data_struct(DataStruct,msSpecs)
%% creating relevant struct for all ch
SpeedOfLight = 299792458;
wavelength = SpeedOfLight/msSpecs.freq0;
Baselines = sqrt(msSpecs.uvw(:,1).^2+msSpecs.uvw(:,2).^2);
flag = (~DataStruct.flag(:)).*(DataStruct.weights(:)>0).*(abs(DataStruct.y_I(:))>0) .*(Baselines(:)>0);
clear Baselines;
% Time related
time = msSpecs.ObsTime(:) ; msSpecs.time=[];
time = time(flag>0);
gSnapNum  = nnz(unique(time));
snapshots = unique(time);
[~,timeSnapIdx] = ismember(snapshots,time);
timeSnapIdx=[timeSnapIdx;length(time)];
% Natural weights
nW = DataStruct.weights;
nW = sqrt(nW(:));
nW = double(nW(flag>0));
% Apply weights to data
y  = double((DataStruct.y_I(:)));
y  = y(flag>0).*nW;
% Normalise uvw coordinates in units of the wavelength
u = msSpecs.uvw(:,1)./wavelength;
v = msSpecs.uvw(:,2)./wavelength;
w = msSpecs.uvw(:,3)./wavelength;
u = u(flag>0);
v = v(flag>0);
w = w(flag>0);

% Antenna related
ant1     = msSpecs.ant1(:) ;msSpecs.ant1=[];
ant2     = msSpecs.ant2(:) ;msSpecs.ant2=[];
ant1 = ant1(flag>0);
ant2 = ant2(flag>0);
ants21   = [ant2(:) ant1(:)];
gAntNum  = max(max(ants21));
AntsPerSnapNum = zeros(gSnapNum,1);
for timeSnap= 1: gSnapNum
    AntsPerSnapNum(timeSnap)=nnz(unique(ants21(timeSnapIdx(timeSnap):timeSnapIdx(timeSnap+1)-1,:)));
end

% Update Flag col
flag = flag(flag>0) ;

%% Building struct. and re-ordering data

% snap_gen.timeInstants = unique(time);
snap_gen.active_ants_nber = (nnz(unique(ants21(:))));
snap_gen.AntsPerSnapNum = AntsPerSnapNum;

% Init
snap_gen.time = cell(gSnapNum,1);
snap_gen.ant1 = cell(gSnapNum,1);
snap_gen.ant2 = cell(gSnapNum,1);
snap_gen.u = cell(gSnapNum,1);
snap_gen.v = cell(gSnapNum,1);
snap_gen.w = cell(gSnapNum,1);
%
data_spec.y    = cell(gSnapNum,1);
data_spec.flag = cell(gSnapNum,1);
data_spec.nW   = cell(gSnapNum,1);
%
snap_gen.start = zeros(gSnapNum,1);
snap_gen.end    = zeros(gSnapNum,1);
snap_gen.sizes  = zeros(gSnapNum,1);
snap_gen.involved_ants_idx = cell(gSnapNum,1);
snap_gen.active_ants = cell(gSnapNum,1);
snap_gen.active_ants_nber =  zeros(gSnapNum,1);
snap_gen.order1 = cell(gSnapNum,1);
snap_gen.order2 = cell(gSnapNum,1);
%
snap_gen.nonActiveAnts =[];
%
fprintf('\nRe-ordering data & preparing structures ..')
for timeSnap = 1:gSnapNum
    idStart = timeSnapIdx(timeSnap);
    if timeSnap <gSnapNum
        idEnd   = timeSnapIdx(timeSnap+1)-1;
    else, idEnd   = timeSnapIdx(end);
    end
    snap_gen.start(timeSnap) = idStart;
    snap_gen.end(timeSnap)   = idEnd;
    flagS = flag(idStart: idEnd);
    if (nnz(flagS))
        unaS = u(idStart: idEnd);
        vnaS = v(idStart: idEnd);
        wnaS = w(idStart: idEnd);
        na1S = ants21(idStart: idEnd,1) ;
        na2S = ants21(idStart: idEnd,2) ;
        timeS = time(idStart: idEnd);
        dataS = y(idStart: idEnd);
        nWS = nW(idStart: idEnd);
        
        order1 =[]; order2 =[];
        
        [~,order1] = sort(na1S,'descend');
        na1S = na1S(order1);
        na2S = na2S(order1);
        unaS = unaS(order1);
        vnaS = vnaS(order1);
        wnaS = wnaS(order1);
        timeS = timeS(order1);
        dataS = dataS(order1);
        nWS   = nWS(order1);
        flagS = flagS(order1);
        %
        order2 =[];
        for step =max(unique(na1S)):-1:min(unique(na1S))
            id = find(~(na1S-step));
            [~, v_sort] = sort(na2S(id),'descend');
            na2S(id)   = na2S(id(v_sort));
            na1S(id)   = na1S(id(v_sort));
            timeS(id)  = timeS(id(v_sort));
            unaS(id)   = unaS(id(v_sort));
            vnaS(id)   = vnaS(id(v_sort));
            wnaS(id)   = wnaS(id(v_sort));
            dataS(id)  = dataS(id(v_sort));
            flagS(id)  = flagS(id(v_sort));
            nWS(id)    = nWS(id(v_sort));
            order2     = [order2; id(v_sort)]; %#ok<AGROW>
        end
        
        snap_gen.time{timeSnap} = timeS;
        snap_gen.ant1{timeSnap} = na1S;
        snap_gen.ant2{timeSnap} = na2S;
        snap_gen.u{timeSnap} = unaS;
        snap_gen.v{timeSnap} = vnaS;
        snap_gen.w{timeSnap} = wnaS;
        snap_gen.sizes(timeSnap)= nnz(na1S);
        data_spec.y{timeSnap}  = dataS;
        data_spec.flag{timeSnap}  = flagS;
        data_spec.nW{timeSnap}  = nWS;
        [AntsSnapInvolved,~,~]=unique(nonzeros(ants21(idStart: idEnd,:)));
        tmp = zeros(gAntNum,1);
        for j=1:length(AntsSnapInvolved)
            tmp(AntsSnapInvolved(j))=j;      
        end
        snap_gen.involved_ants_idx{timeSnap} = tmp;
        snap_gen.active_ants{timeSnap} = AntsSnapInvolved;
        snap_gen.active_ants_nber(timeSnap)=length(AntsSnapInvolved);
        %
        snap_gen.order1{timeSnap}= order1;
        snap_gen.order2{timeSnap}= order2; 
        
    else,    snap_gen.nonActiveAnts=[snap_gen.nonActiveAnts;timeSnap];
    end
end
%
data_spec.nW   = cell2mat(data_spec.nW );
gVars.gAntNum  = gAntNum;
gVars.gSnapNum = gSnapNum;
gVars.MeasNum  = numel(flag);
end

