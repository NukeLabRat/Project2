uCentral=interp2(uXlocations,uYlocations,u(:,:,k),pXlocations,pYlocations);
vCentral=interp2(vXlocations,vYlocations,v(:,:,k),pXlocations,pYlocations); 
yLines=25:4:45;
figure
Plot1(1)=plot(uCentral(:,floor(Nodes/2)),pYlocations(:,floor(Nodes/2)),'bo','DisplayName','Gould','MarkerFaceColor','b','MarkerSize',12);
Bench=csvread('KimRe400Profile31.csv',1);
hold on
Plot1(2)=plot(Bench(:,1),Bench(:,2),'rsq','MarkerFaceColor','r','MarkerSize',13,'DisplayName','Ghia (1982)');
axes3=gca;
set(axes3,'FontName','Times New Roman','FontSize',28,'Layer','top','LineWidth',3);
xlabel('u-velocity')
ylabel('y')
title('Velocity Profiles - Re 1000')
grid on
axes3.Color='none';
ylim([0 1])