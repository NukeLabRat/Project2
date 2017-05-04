close all
figure;
uCentral(isnan(uCentral)==1)=0;
vCentral(isnan(vCentral)==1)=0;
[phi, psi]=flowfun(uCentral,vCentral);
contour(pXlocations,pYlocations,phi,20,'LineWidth',3);
colorbar
axes1=gca;
set(axes1,'FontName','Times New Roman','FontSize',28,'LineWidth',3);
title('\phi')
xlabel('x')
ylabel('y')
axes1.Color='none';
figure
contour(pXlocations,pYlocations,psi,20,'LineWidth',3);
colorbar
axes4=gca;
set(axes4,'FontName','Times New Roman','FontSize',28,'LineWidth',3);
xlabel('x')
ylabel('y')
title('\psi')
figure
axes4.Color='none';
pcolor(pXlocations(2:end-1,2:end-1),pYlocations(2:end-1,2:end-1),uCentral(2:end-1,2:end-1));

axes3=gca;

% % colorbar
% yy=1:Nodes;xx=1:Nodes;
% % Horizontal grid
% 
% for k = 1:length(yy)
% line([xx(1) xx(end)], [yy(k) yy(k)],'Color',[0 0 0],'LineWidth',1)
% end
% % Vertical grid
% for k = 1:length(xx)
% line([xx(k) xx(k)], [yy(1) yy(end)],'Color',[0 0 0],'LineWidth',1)
% end

% 
% for k = 1:length(y)
% line([x(1) x(end)], [y(k) y(k)],'Color',[0 0 0],'LineWidth',1)
% end
% % Vertical grid
% for k = 1:length(x)
% line([x(k) x(k)], [y(1) y(end)],'Color',[0 0 0],'LineWidth',1)
% end
% xlim([0 Nodes-1])
% ylim([0 Nodes-1])
% NodeInterp=
NodeLabelsNum=floor(linspace(1,50,10));
NodeLabels=num2cell(NodeLabelsNum);
for i=1:length(NodeLabels)
    Ninterp=interp1(NodeLabelsNum,pXlocations(1,NodeLabelsNum));
NodeLabels{i}=num2str(NodeLabels{i});
end
% set(axes3,'FontName','Times New Roman','FontSize',25,'Layer','top','LineWidth',3,'XTick',NodeLabelsNum,...
%     'XTickLabel',NodeLabels,'YTick',NodeLabelsNum,'YTickLabel',NodeLabels);
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
