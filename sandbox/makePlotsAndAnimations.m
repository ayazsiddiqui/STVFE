%% plot data
% % % linewidths
lwd = 1;
% % % set find plot limits
lB = min([lowerBound(:);windSpeedOut(:)]);
uB = max([upperBound(:);windSpeedOut(:)]);
plotRes = 1;

figure(1)
set(gcf,'Units','normalized','position',[1 0.0889 0.6667 0.8398]);
F = struct('cdata',uint8(zeros(840,1680,3)),'colormap',[]);

for ii = 1:noTimeSteps
    
    if ii == 1
        hold on
        grid on
        xlabel('Wind speed (m/s)');
        ylabel('Altitude (m)');
        %         xlim([lB-mod(lB,plotRes),uB-mod(uB,plotRes)+plotRes])
        xlim(meanFlow+[-4 4])
        ylim([hMin hMax]);
    else
        delete(findall(gcf,'type','annotation'));
        h = findall(gca,'type','line','linestyle','-','-or',...
            'linestyle','--','-or','linestyle','-x','-or',...
            'linestyle','-.', '-or','color','m');
        delete(h);
        
    end
    %     end
    
    % % plot true wind
    %     for jj = 1:predHorz
    %         subplot(2,2,jj)
    plTrueWind = plot(windSpeedOut(:,ii),heights,'k','linewidth',lwd);
    % % plot measured wind value
    plfVals = plot(fValAtPt(ii,:),pointsVisited(:,ii),'mo',...
        'markerfacecolor','m','linewidth',lwd);
    % % plot GPKF mean and bounds
    plPredMean = plot(predMean(:,ii),xPredict,'-x','linewidth',lwd,...
        'color',1/255*[228,26,28]);
    plLowerBds = plot(lowerBound(:,ii),xPredict,'--','linewidth',lwd,...
        'color',1/255*[254,178,76]);
    plUpperBds = plot(upperBound(:,ii),xPredict,'--','linewidth',lwd,...
        'color',1/255*[254,178,76]);
    % % legend
    legend([plTrueWind,plPredMean,plLowerBds],...
        'True func','GPKF $\mu$','GPKF bounds');
    % % title
    txt1 = sprintf('$l_{t} / \\tau$ = %0.2f,',timeScale/timeStep);
    txt = sprintf(' Time = %0.2f min',tVec(ii));
    txt = strcat(txt1,txt);
    title(txt);
    %     end
    ff = getframe(gcf);
    F(ii).cdata = ff.cdata;
    
end

%% other plots
figure(2)
x = gcf;
lwd2 = 1;
set(gcf,'position',x.Position.*[1 0 1 1])
hold on
grid on
pFvalBF = stairs(0:mpcCount-2,optFvalBF,'--o','linewidth',lwd2,...
    'color',1/255*[228,26,28]);
pFvalFmin = stairs(1:mpcCount-1,optFvalFmin,'--o','linewidth',lwd2,...
    'color',1/255*[254,178,76]);
xlabel('MPC run count')
ylabel('Objective function value')
legend([pFvalBF,pFvalFmin],{'Brute force','SQP'})
set(findobj('-property','FontSize'),'FontSize',12)

%% video
% % % video setting
video = VideoWriter(strcat(fName,'_windProf'),'Motion JPEG AVI');
% % video = VideoWriter('vid_Test1','MPEG-4');
video.FrameRate = 5;
set(gca,'nextplot','replacechildren');

open(video)
for ii = 1:length(F)
    writeVideo(video, F(ii));
end
close(video)

%%
figure(3)
set(gcf,'position',x.Position.*[1 0 1 1])
F = struct('cdata',uint8(zeros(840,1680,3)),'colormap',[]);

for ii = 1:size(optStateTrajBF,1)
    
    if ii == 1
        hold on
        grid on
        xlabel('Prediction horizon');
        ylabel('Altitude (m)');
        %         xlim([lB-mod(lB,plotRes),uB-mod(uB,plotRes)+plotRes])
        xlim([0 predHorz+1])
        ylim([hMin hMax]);
    else
        delete(findall(gcf,'type','annotation'));
        h = findall(gca,'type','line','linestyle','-','-or',...
            'linestyle','--','-or','linestyle','-x','-or',...
            'linestyle','-.', '-or','color','m');
        delete(h);
        
    end
    %     end
    
    % % plot true wind
    %     for jj = 1:predHorz
    %         subplot(2,2,jj)
    plBF = stairs(0:predHorz,[initCon(ii) optStateTrajBF(ii,:)],'r--o',...
        'linewidth',lwd);
    % % plot measured wind value
    plFmin = stairs(0:predHorz,[initCon(ii) optStateTrajFmin(ii,:)],'b--o',...
        'linewidth',lwd);
    % % legend
    legend([plBF,plFmin],...
        'Brute force','SQP');
    % % title
    txt1 = sprintf('MPC count:');
    txt = sprintf(' %d',ii);
    txt = strcat(txt1,txt);
    title(txt);
    %     end
    ff = getframe(gcf);
    F(ii).cdata = ff.cdata;
    
end

%% video
% % % video setting
video = VideoWriter(strcat(fName,'_optSeq'),'Motion JPEG AVI');
% % video = VideoWriter('vid_Test1','MPEG-4');
video.FrameRate = 5;
set(gca,'nextplot','replacechildren');

open(video)
for ii = 1:length(F)
    writeVideo(video, F(ii));
end
close(video)