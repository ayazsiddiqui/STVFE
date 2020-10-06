function val = parseLogsout(inVal,varargin)

pp = inputParser;
addParameter(pp,'resampleTime',0.25,@isnumeric);

parse(pp,varargin{:});

%
val = struct;

% store in timeseries structure
for ii = 1:inVal.logsout.numElements
    % processed signal names
    pNames = erase(inVal.logsout{ii}.Name,["<",">"," "]);
    val.(pNames) = resample(inVal.logsout{ii}.Values,...
        0:pp.Results.resampleTime:inVal.tout(end));

end