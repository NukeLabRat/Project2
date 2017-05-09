uCentral=interp2(uXlocations,uYlocations,u(:,:,k),pXlocations,pYlocations);
vCentral=interp2(vXlocations,vYlocations,v(:,:,k),pXlocations,pYlocations); 
yLines=25:4:45;
figure
plot(vCentral(:,yLines),pYlocations(:,yLines))
Bench=csvread('Benchmark.csv',1);
hold on
plot(Bench(:,1),Bench(:,2),'bsq','MarkerFaceColor','b','MarkerSize',10)