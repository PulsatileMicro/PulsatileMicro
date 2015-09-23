%% Get .his data file according to given parameters
function []=GetHisFile(VesType, VesParam, BCType, Bifur, BCVal, dt, Nstep)

for i=2:2
    switch i
        case 1
            HisPos='start';
        case 2
            HisPos='mid';
        case 3
            HisPos='end';
    end
    fileName=GenInput4Network(VesType, VesParam, HisPos, BCType, Bifur, BCVal, dt, Nstep);
    system(['./oneDbio_visc ' fileName]);
end

