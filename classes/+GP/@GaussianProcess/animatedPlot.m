function val = animatedPlot(obj,tsDataStruct)

% flow time series object
flowTs = tsDataStruct.FlowData;
% altitude time series object
altTs = tsDataStruct.altData;
% resample
% time values
tVals = flowTs.Time;
% number of time steps
nTs = numel(tVals);

% local variables
flowVals = flowTs.Data;
altVals  = altTs.Data;

% create axis object
axisObj = axes;
axisObj.YLim = [altVals(1) altVals(end)];

% create plot object
plotObj = plot(altVals,flowVals(:,:,1));



for ii = 2:nTs
    
    flowVals.YData = flowVals(:,:,ii);
    
    
    
    
end
end

