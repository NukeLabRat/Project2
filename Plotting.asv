close all
figure;
k=1;
uCentral=interp2(uXlocations,uYlocations,u(:,:,k),pXlocations,pYlocations);
vCentral=interp2(vXlocations,vYlocations,v(:,:,k),pXlocations,pYlocations);
uCentral(isnan(uCentral)==1)=0;
vCentral(isnan(vCentral)==1)=0;
[phi, psi]=flowfun(uCentral,vCentral);


pcolor(pXlocations(2:end-1,2:end-1),pYlocations(2:end-1,2:end-1),(uCentralEnd(2:end-1,2:end-1).^2+vCentralEnd(2:end-1,2:end-1).^2).^.5);
SpeedPlot = gca;
box(SpeedPlot,'on');
set(SpeedPlot,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
colorbar
title('Flow Speed')
SpeedPlot.TickLength=[0 0];
SpeedPlot.YDir='normal';


figure
contour(pXlocations,pYlocations,psi,40,'LineWidth',3,'LineColor','b');
[StreamX, StreamY]=meshgrid([dx .1 .15 .2 .85 .88 .9 .95],dy*ones(1,10));
XY=stream2(pXlocations,pYlocations,uCentral,vCentral,StreamX,StreamY);
StreamA=streamline(XY)
% axes5=gca;
% set(axes5,'FontName','Times New Roman','FontSize',28,'LineWidth',3)
title('Streamlines')

% colorbar
axes4=gca;
set(axes4,'FontName','Times New Roman','FontSize',28,'LineWidth',3);
xlabel('x')
ylabel('y')
for i=1:length(StreamA)
StreamA(i).LineWidth=3;
end
axes4.Color='none';

figure
pcolor(pXlocations(2:end-1,2:end-1),pYlocations(2:end-1,2:end-1),uCentral(2:end-1,2:end-1));
axes3=gca;
set(axes3,'FontName','Times New Roman','FontSize',28,'Layer','top','LineWidth',3);
xlabel('x')
ylabel('y')
title('u-velocity')
axes3.TickLength=[0 0];
colorbar


figure;
pcolor(pXlocations(2:end-1,2:end-1),pYlocations(2:end-1,2:end-1),vCentral(2:end-1,2:end-1));
axes2 = gca;
box(axes2,'on');
set(axes2,'FontName','Times New Roman','FontSize',28,'LineWidth',3)
colorbar
title('v-velocity')
axes2.YDir='normal';
y=1:Nodes;x=1:Nodes;
% Horizontal grid
% for k = 1:length(y)
% line([x(1) x(end)], [y(k) y(k)],'Color',[0 0 0],'LineWidth',1)
% end
% % Vertical grid
% for k = 1:length(x)
% line([x(k) x(k)], [y(1) y(end)],'Color',[0 0 0],'LineWidth',1)
% end
% xlim([1 Nodes])
% ylim([1 Nodes])
% NodeLabelsNum=1:6:Nodes;
% NodeLabels=num2cell(1:6:Nodes);
% for i=1:length(NodeLabels)
% NodeLabels{i}=num2str(NodeLabels{i});
% end
% set(axes2,'FontName','Times New Roman','FontSize',25,'Layer','top','LineWidth',3,'XTick',NodeLabelsNum,...
%     'XTickLabel',NodeLabels,'YTick',NodeLabelsNum,'YTickLabel',NodeLabels);
set(axes2,'FontName','Times New Roman','FontSize',28,'Layer','top','LineWidth',3);
xlabel('x')
ylabel('y')
axes2.TickLength=[0 0];

% figure
% [StreamX, StreamY]=meshgrid([dx .05 .1 .15 .2 .22 .85 .88 .9 .95],dy*ones(1,10));
% XY=stream2(pXlocations,pYlocations,uCentral,vCentral,StreamX,StreamY);
% streamline(XY)
% axes5=gca;
% set(axes2,'FontName','Times New Roman','FontSize',28,'LineWidth',3)
% title('Streamlines')

