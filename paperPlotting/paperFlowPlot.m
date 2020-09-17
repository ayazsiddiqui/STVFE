function axisObj = paperFlowPlot(synFlow,altitudes,endTime)

%% process data
% plot struct
P.plotTimeStep = 10;
P.endTime      = endTime;
P.xlimRes      = 2;
cols = [228,26,28
    55,126,184
    77,175,74
    152,78,163
    255,127,0
    166,86,40
    231,41,138]/255;
P.linProp = {...
    cols(1,:),'-';
    cols(2,:),':';
    cols(3,:),'-^';
    cols(4,:),'--';
    cols(5,:),'-o';
    cols(6,:),'-s';
    cols(7,:),'-d'};


tVec = synFlow.Time(1):P.plotTimeStep*60:P.endTime*60;
flowTs = resample(synFlow,tVec);
lb = min(flowTs.Data,[],'all');
ub = max(flowTs.Data,[],'all');

P.xLoLim = floor(lb) - mod(floor(lb),P.xlimRes);
P.xHiLim = ceil(ub) + mod(ceil(ub),P.xlimRes);

% create axis object
axisObj = axes;
% set axis properties
grid(axisObj,'on');
set(gca,'GridLineStyle',':')
hold(axisObj,'on');
xlabel(axisObj,'\textbf{Flow speed [m/s]}','fontweight','bold');
ylabel(axisObj,'\textbf{Altitude [m]}','fontweight','bold');
set(axisObj,'FontSize',12);
axisObj.YLim = [altitudes(1) altitudes(end)];
% set value
axisObj.XLim = [P.xLoLim P.xHiLim];
set(gcf,'InnerPosition',1*[-0 -0 560 420])

%% plot the data
fPlots = gobjects;
for ii = 1:length(tVec)
    legStr{ii} = strcat(num2str(tVec(ii)/60),' min');
    fPlots(ii) = plot(flowTs.Data(:,:,ii),altitudes...
        ,P.linProp{ii,2}...
        ,'color',P.linProp{ii,1}...
        ,'MarkerFaceColor',P.linProp{ii,1}...
        ,'MarkerSize',5 ...
        ,'linewidth',1);
    hold on
end

legend(fPlots(:),legStr{:},'location','southeast')

end