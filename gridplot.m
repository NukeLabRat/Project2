for i=1:7
for j=1:6
Vtext{j,i}=['$$v_{' num2str(i-1) ',' num2str(j) '\frac{1}{2}}$$'];
end
end

for i=1:6
for j=1:7
Utext{j,i}=['$$u_{' num2str(i) '\frac{1}{2},' num2str(j-1) '}$$'];
end
end

for i=1:7
for j=1:7
Ptext{j,i}=['$$P_{' num2str(i-1) ',' num2str(j-1) '}$$'];
end
end
figure
hold on

plot(uXlocations(1),uYlocations(1),'r.','MarkerSize',15,'DisplayName','u-velocity')
plot(vXlocations(1),vYlocations(1),'b.','MarkerSize',15)
plot(pXlocations(1),pYlocations(1),'k.','MarkerSize',15)
Axes1=gca;
% Legend1=legend(Axes1,'show');

set(Axes1, 'color', 'none');
Axes1.Visible='off';
% set(Legend1, 'color', 'none');
box off
set(Axes1,'FontName','Times New Roman','FontSize',18,'LineWidth',3)
line([0 1], [1 1],'Color',[81/255 40/255 136/255],'LineWidth',2)
line([1 1], [0 1],'Color',[81/255 40/255 136/255],'LineWidth',2)
line([0 1], [0 0],'Color',[81/255 40/255 136/255],'LineWidth',2)
line([0 0], [0 1],'Color',[81/255 40/255 136/255],'LineWidth',2)
plot(uXlocations,uYlocations,'r.','MarkerSize',15)
plot(vXlocations,vYlocations,'b.','MarkerSize',15)
plot(pXlocations,pYlocations,'k.','MarkerSize',15)


for i=1:6
for j=1:7
    
text(uXlocations(j,i)-.035,uYlocations(j,i)+.035,Utext{j,i},'interpreter','latex','FontSize',17,'Color','r')

end
end

for i=1:7
for j=1:6
text(vXlocations(j,i)-.035,vYlocations(j,i)+.035,Vtext{j,i},'interpreter','latex','FontSize',17,'Color','b')
end
end

for i=1:7
for j=1:7
text(pXlocations(j,i)-.035,pYlocations(j,i)+.035,Ptext{j,i},'interpreter','latex','FontSize',17,'Color','k')
end
end


xlim([-1.5*dx Nodes*dx+.5*dx])
ylim([-1.5*dy Nodes*dy+.5*dy])

set(Axes1,'FontName','Times New Roman','FontSize',25,'Layer','top','LineWidth',3,'XTick',[],...
    'XTickLabel',[],'YTick',[],'YTickLabel',[]);

Axes1.TickLength=[0 0];
