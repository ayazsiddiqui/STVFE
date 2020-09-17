clear
clc


altitudes = 100:100:500;
times     = 1:4;

[x,y] = meshgrid(times,altitudes);

Mk = [100;400;200;300];

for ii = 1:length(altitudes)
yT{ii} = sprintf('$z_{%d}$',ii);
end

for ii = 1:length(times)+1
xT{ii} = sprintf('$t_{%d}$',ii);
end

axisObj = axes;
% set axis properties
grid(axisObj,'on');
set(gca,'GridLineStyle',':')
hold(axisObj,'on');
xlabel(axisObj,'${t}$','fontweight','bold');
ylabel(axisObj,'$\mathbf{X}$','fontweight','bold');

mSize = 10;
pX = plot(x(:),y(:),'rx','MarkerSize',mSize,'linewidth',1.2);
pM = plot(times,Mk,'ko','MarkerSize',mSize/2,'MarkerFaceColor','k','linewidth',1.2);
set(gca,'ytick',altitudes,'yticklabel',yT)
set(gca,'xtick',[times times(end)+1],'xticklabel',xT)
pP = plot(times(end)+1,250,'mh','MarkerSize',mSize,...
    'MarkerFaceColor','m','MarkerEdgeColor','k','linewidth',1.5);
plot([0 times(end)+1],250*[1 1],'k-.');
plot((times(end)+1)*[1 1],250*[0 1],'k-.');

ylim([0 575])
xlim([0 5.5])
set(gcf,'InnerPosition',1*[-0 -0 560 420])
legend([pX;pM;pP],{'$\mathbf{X_{\textrm{M}}}$','measurements','prediction'},...
    'location','bestoutside','orientation','horizontal')
set(axisObj,'FontSize',15);

 % determine position of the axes

keyboard
% determine startpoint and endpoint for the arrows 
axp = get(gca,'Position');
xs=axp(1);
xe=axp(1)+axp(3);
ys=axp(2);
ye=axp(2)+axp(4);
% make the arrows
annotation('arrow', [xs xe],[ys ys]);
annotation('arrow', [xs xs],[ys ye]);

%% export file
saveFile = input('Save file? Options: Enter y or n\n','s');
if strcmpi(saveFile,'y')
    filName = strcat('kfgpDemo_',strrep(datestr(datetime),':','-'));
    savefig(filName);
    exportgraphics(gcf,[filName,'.png'],'Resolution',600)
end



