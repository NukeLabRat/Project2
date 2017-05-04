close all
figure;
[phi, psi]=flowfun(uCentral(2:end-1,2:end-1),vCentral(2:end-1,2:end-1));
contour(pXlocations(2:end-1,2:end-1),pYlocations(2:end-1,2:end-1),psi,20,'LineWidth',3,'LineColor',[0 0 1]);
imagesc(uCentral(2:end-1,:));
axes3 = gca;
box(axes3,'on');
set(axes3,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
title('\psi')
axes3.YDir='normal';
% colorbar
% y=1:Nodes;x=1:Nodes;
% Horizontal grid
for k = 1:length(y)
line([x(1) x(end)], [y(k) y(k)],'Color',[0 0 0],'LineWidth',1)
end
% Vertical grid
for k = 1:length(x)
line([x(k) x(k)], [y(1) y(end)],'Color',[0 0 0],'LineWidth',1)
end
% xlim([0 Nodes-1])
% ylim([0 Nodes-1])
% NodeInterp=
NodeLabelsNum=1:6:Nodes;
NodeLabels=num2cell(NodeLabelsNum);
for i=1:length(NodeLabels)
NodeLabels{i}=num2str(NodeLabels{i});
end
set(axes3,'FontName','Times New Roman','FontSize',25,'Layer','top','LineWidth',3);
xlabel('x')
ylabel('y')
axes3.TickLength=[0 0];

figure;
imagesc(vCentral(:,2:end-1));
axes2 = gca;
box(axes2,'on');
set(axes2,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
colorbar
title('v')
axes2.YDir='normal';
y=1:Nodes;x=1:Nodes;
% Horizontal grid
for k = 1:length(y)
line([x(1) x(end)], [y(k) y(k)],'Color',[0 0 0],'LineWidth',1)
end
% Vertical grid
for k = 1:length(x)
line([x(k) x(k)], [y(1) y(end)],'Color',[0 0 0],'LineWidth',1)
end
xlim([1 Nodes])
ylim([1 Nodes])
NodeLabelsNum=1:6:Nodes;
NodeLabels=num2cell(1:6:Nodes);
for i=1:length(NodeLabels)
NodeLabels{i}=num2str(NodeLabels{i});
end
set(axes2,'FontName','Times New Roman','FontSize',25,'Layer','top','LineWidth',3,'XTick',NodeLabelsNum,...
    'XTickLabel',NodeLabels,'YTick',NodeLabelsNum,'YTickLabel',NodeLabels);
xlabel('x')
ylabel('y')
axes2.TickLength=[0 0];
