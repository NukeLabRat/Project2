%% CFD PROJECT 2
% Dan Gould

close all
clear
clc
StartTime=tic;
Re=800;
dt=.004;
Nodes=50;
AllowedError=5E-8;
PoissonIn.SOR=.95;
PoissonIn.Error=5E-10;

%% Geometry -
L = 1; %m, y-dir
W = 1; %m, x-dir

%%
dx = 1/(Nodes-1); %m
dy = 1/(Nodes-1); %m
Beta = dx/dy;

x = 0:dx:W; %vector of node locations; x-dir
y = 0:dy:L; %vector of node locations; y-dir
xPlusHalf = x+.5*dx;
yPlusHalf = y+.5*dy;
xEnd = length(x)+1; %length of location vectors
yEnd = length(y)+1;

[X, Y] = meshgrid(x,y);
[pXlocations, pYlocations] = meshgrid([-.5*dx xPlusHalf],[-.5*dy yPlusHalf]);
[uXlocations, uYlocations] = meshgrid(x,[-.5*dy yPlusHalf]);
[vXlocations, vYlocations] = meshgrid([-.5*dx xPlusHalf],y);
uSize=size(uXlocations);
vSize=size(vXlocations);
pSize=size(pXlocations);

P = zeros(yEnd,xEnd); %initialize matrix with zeros
P0=P;
NodesP = 1:numel(P);
NodesX = 1:numel(uXlocations);
NodesY = 1:numel(uXlocations);
%% Stability
Stab1=.25*dt*Re
Stab2=dt/(Re*dx^2)

if Stab1>.95 || Stab2>.25
    error('Unstable')
end

%% Apply BC
Pinitial = 0;
BoundaryNodesP=MatEdges(P);
BoundaryNodesX=MatEdges(uXlocations);
BoundaryNodesY=MatEdges(vYlocations);

%% Center Nodes
CenterNodesP = NodesP(~ismember(NodesP,BoundaryNodesP));
CenterNodesX = NodesX(~ismember(NodesX,BoundaryNodesX));
CenterNodesY = NodesY(~ismember(NodesY,BoundaryNodesY));

IndexX=zeros(uSize(1),uSize(2));
IsCenterX=zeros(uSize(1),uSize(2));
for i = 1:uSize(2) %This loop determines whether each node is central node or boundary node.
    for j = 1:uSize(1)
        IndexX(j,i) = sub2ind([uSize(1) uSize(2)],j,i);
        if sum(ismember(CenterNodesX,IndexX(j,i)))==1
            IsCenterX(j,i) = true;
        else
            IsCenterX(j,i) = false;
        end
    end
end
IndexY=zeros(vSize(1),vSize(2));
IsCenterY=zeros(vSize(1),vSize(2));
for i = 1:vSize(2) %This loop determines whether each node is central node or boundary node.
    for j = 1:vSize(1)
        IndexY(j,i) = sub2ind([vSize(1) vSize(2)],j,i);
        if sum(ismember(CenterNodesY,IndexY(j,i)))==1
            IsCenterY(j,i) = true;
        else
            IsCenterY(j,i) = false;
        end
    end
end

IndexP=zeros(yEnd,xEnd);
IsCenterP=zeros(yEnd,xEnd);
for i = 1:xEnd %This loop determines whether each node is central node or boundary node.
    for j = 1:yEnd
        IndexP(j,i) = sub2ind([yEnd xEnd],j,i);
        if sum(ismember(CenterNodesP,IndexP(j,i)))==1
            IsCenterP(j,i) = true;
        else
            IsCenterP(j,i) = false;
        end
    end
end
IsCenterP=logical(IsCenterP);
TopWallP=false(pSize);
BottomWallP=false(pSize);
RightWallP=false(pSize);
LeftWallP=false(pSize);
    for i = 1:xEnd
        for j = 1:yEnd
            if IsCenterP(j,i)==false %checks if node is central node
                if j==1
                    BottomWallP(j,i)=true;
%                     Pold(j,i)=Pold(j+1,i);
                end
                if j==yEnd
                    TopWallP(j,i)=true;
%                     Pold(j,i)=Pold(j-1,i);
                end
                if i==1
                    LeftWallP(j,i)=true;
%                     Pold(j,i)=Pold(j,i+1);
                end
                if i==xEnd
                    RightWallP(j,i)=true;
%                     Pold(j,i)=Pold(j,i-1);
                end
            end   
        end 
    end 
    PoissonIn.Jminus=[];
    PoissonIn.Jplus=[];
    PoissonIn.Iminus=[];
    PoissonIn.Iplus=[];
    PoissonIn.Center=[];
    
     for i = 1:xEnd
        for j = 1:yEnd
            if IsCenterP(j,i)==true
  PoissonIn.Center(end+1)=sub2ind([yEnd xEnd],j,i);              
  PoissonIn.Jminus(end+1)=sub2ind([yEnd xEnd],j-1,i);
  PoissonIn.Jplus(end+1)=sub2ind([yEnd xEnd],j+1,i);
  PoissonIn.Iminus(end+1)=sub2ind([yEnd xEnd],j,i-1);
  PoissonIn.Iplus(end+1)=sub2ind([yEnd xEnd],j,i+1);
            end
        end 
     end 
     
PoissonIn.IsCenterP=IsCenterP;
PoissonIn.TopWallP=false(pSize);
PoissonIn.BottomWallP=false(pSize);
PoissonIn.RightWallP=false(pSize);
PoissonIn.LeftWallP=false(pSize);
PoissonIn.TopWallPmirror=false(pSize);
PoissonIn.BottomWallPmirror=false(pSize);
PoissonIn.RightWallPmirror=false(pSize);
PoissonIn.LeftWallPmirror=false(pSize);

    for i = 1:xEnd
        for j = 1:yEnd
            if IsCenterP(j,i)==false %checks if node is central node
                if j==1
                    PoissonIn.BottomWallP(j,i)=true;
                    PoissonIn.BottomWallPmirror(j+1,i)=true;
                end
                if j==yEnd
                    PoissonIn.TopWallP(j,i)=true;
                    PoissonIn.TopWallPmirror(j-1,i)=true;
                end
                if i==1
                    PoissonIn.LeftWallP(j,i)=true;
                    PoissonIn.LeftWallPmirror(j,i+1)=true;
                end
                if i==xEnd
                    PoissonIn.RightWallP(j,i)=true;
                    PoissonIn.RightWallPmirror(j,i-1)=true;
                end
            end   
        end 
    end %Create 
    

PoissonIn.dx=dx;
PoissonIn.dy=dy;
PoissonIn.xSize=pSize(2);
PoissonIn.ySize=pSize(1);
TimeLevals=2;
u=zeros(uSize(1),uSize(2),TimeLevals);
v=zeros(vSize(1),vSize(2),TimeLevals);
P=zeros(pSize(1),pSize(2),TimeLevals);
ConstantMat=zeros(pSize(1),pSize(2));
Ustar(:,:,1)=ones(uSize(1),uSize(2));
Vstar(:,:,1)=ones(vSize(1),vSize(2));
dxUstar(:,:,1)=ones(uSize(1),uSize(2));
dyVstar(:,:,1)=ones(vSize(1),vSize(2));
duStarCentral(:,:,1)=ones(pSize(1),pSize(2));
dvStarCentral(:,:,1)=ones(pSize(1),pSize(2));
SOR=1;

for i = 1:xEnd-1 %Assign velocity BC for initial Step
    for j = 1:yEnd
        if IsCenterX(j,i)==true %checks if node is central node
        else %For Boundary Nodes
            if j==1
                u(j,i,2)=-u(2,i,2)/2;
            elseif j==uSize(1)
                u(j,i,2)=2-u(uSize(1)-1,i,2);
            else
                u(j,i,2)=0;
            end
        end
    end
end
PressureIterations=zeros(1,1000);
TimeCheck=0;
Error2=1;
MainIterations=1;
while Error2>AllowedError || MainIterations<100
    
    if Error2<.05
        PoissonIn.SOR=1.2;
    else
        PoisonIn.SOR=.7;
    end
    
        
    
    u(:,:,1)=u(:,:,2);
    v(:,:,1)=v(:,:,2);
    P(:,:,1)=P(:,:,2);
    k=1;
    unow=u(:,:,1);
    vnow=v(:,:,1);
    uCentral=interp2(uXlocations,uYlocations,u(:,:,k),pXlocations,pYlocations);
    vCentral=interp2(vXlocations,vYlocations,v(:,:,k),pXlocations,pYlocations);
    uvdy=(interp2(uXlocations,uYlocations,u(:,:,k),uXlocations,uYlocations+.5*dy).*interp2(vXlocations,vYlocations,v(:,:,k),uXlocations,uYlocations+.5*dy)...
        -interp2(uXlocations,uYlocations,u(:,:,k),uXlocations,uYlocations-.5*dy).*interp2(vXlocations,vYlocations,v(:,:,k),uXlocations,uYlocations-.5*dy))/dy;
    uvdx=(interp2(uXlocations,uYlocations,u(:,:,k),vXlocations+.5*dx,vYlocations).*interp2(vXlocations,vYlocations,v(:,:,k),vXlocations+.5*dx,vYlocations)...
        -interp2(uXlocations,uYlocations,u(:,:,k),vXlocations-.5*dx,vYlocations).*interp2(vXlocations,vYlocations,v(:,:,k),vXlocations-.5*dx,vYlocations))/dx;
    dxuSquared=(uCentral.^2-uCentral.^2)/dx;
    dyvSquared=(vCentral.^2-vCentral.^2)/dy;
    
    for i = 1:xEnd-1
        for j = 1:yEnd
            if IsCenterX(j,i)==true %checks if node is central node
                du2dx2=(u(j,i+1,k)-2*u(j,i,k)+u(j,i-1,k))./dx^2;
                du2dy2=(u(j+1,i,k)-2*u(j,i,k)+u(j-1,i,k))./dy^2;
                Ustar(j,i) = (-dxuSquared(j,i)-uvdy(j,i)+1./Re.*(du2dx2+du2dy2)).*dt+u(j,i,k);
            else %For Boundary Nodes
%                 Ustar(j,i)=u(j,i,k)+.25*(u(j,i,k)+u(j,i,2))/2;
                Ustar(j,i)=u(j,i,k);
%                 Ustar(j,i)=0;
            end
            
        end
    end
    for i = 1:xEnd
        for j = 1:yEnd-1
            if IsCenterY(j,i)==true %checks if node is central node
                dv2dx2=(v(j,i+1,k)-2*v(j,i,k)+v(j,i-1,k))./dx^2;
                dv2dy2=(v(j-1,i,k)-2*v(j,i,k)+v(j+1,i,k))./dy^2;
                Vstar(j,i) = (-dyvSquared(j,i)-uvdx(j,i)+1./Re.*(dv2dx2+dv2dy2)).*dt+v(j,i,k);
            else %For Boundary Nodes
%                 Vstar(j,i) =v(j,i,1)+.25*(v(j,i,k)+v(j,i,2))/2;
                Vstar(j,i) =v(j,i,1);
%                 Vstar(j,i) =0;
            end
        end
    end
    
%     for i = 1:xEnd-1
%         for j = 1:yEnd
%             if IsCenterX(j,i)==true %checks if node is central node
%                 dxUstar(j,i)=(Ustar(j,i+1)-Ustar(j,i))/dx; %is in cell center
%             else %For Boundary Nodes
%                 dxUstar(j,i)=1;
%             end
%         end
%     end
%     
%     for i = 1:xEnd
%         for j = 1:yEnd-1
%             if IsCenterY(j,i)==true %checks if node is central node
%                 dyVstar(j,i)=(Vstar(j+1,i)-Vstar(j,i))/dy;
%             else %For Boundary Nodes
%                 dyVstar(j,i)=1;
%             end
%         end
%     end
    PoissonIn.ConstantMat=padarray((diff(Ustar(2:end-1,:),1,2)/dx+diff(Vstar(:,2:end-1),1,1)/dy)./dt,[1 1]);
    PoissonIn.P0=P(:,:,k);

    
    [Pressure, PoissonIterations]=PoisonPressureSLOR3(PoissonIn);
    PressureIterations(MainIterations)=PoissonIterations;
    P(:,:,k+1)=Pressure;
    PinterpU=interp2(pXlocations,pYlocations,P(:,:,k+1),uXlocations,uYlocations);
    PinterpV=interp2(pXlocations,pYlocations,P(:,:,k+1),vXlocations,vYlocations);
    for i = 1:xEnd-1
        for j = 1:yEnd
            if IsCenterX(j,i)==true %checks if node is central node
                u(j,i,k+1) = Ustar(j,i)+dt/dx*(PinterpU(j,i+1)-PinterpU(j,i));
            else %For Boundary Nodes
                if j==1
                    u(j,i,k+1)=-u(2,i,k);
                elseif j==uSize(1)
                    u(j,i,k+1)=2-u(uSize(1)-1,i,k);
                else
                    u(j,i,k+1)=0;
                end
                
            end
            
        end
    end
    for i = 1:xEnd
        for j = 1:yEnd-1
            if IsCenterY(j,i)==true %checks if node is central node
                v(j,i,k+1) = Vstar(j,i)+dt/dy*(PinterpV(j+1,i)-PinterpV(j,i));
            else %For Boundary Nodes
                if i==1
                    v(j,i,k+1)=-v(j,2,k);
                elseif i==vSize(2)
                    v(j,i,k+1)=-v(j,vSize(2)-1,k);
                else
                    v(j,i,k+1)=0;
                end
            end
        end
    end
    Error2=norm(u(:,:,2)-u(:,:,1));
    MainIterations=MainIterations+1;
TimeCheck=TimeCheck+1;
if TimeCheck==20;
    clc
    Stab1
    Stab2
    CurrentTime=dt*MainIterations
    Error=Error2
    TimeCheck=0;
end
end
StopTime=toc(StartTime)
PressureIterations(PressureIterations==0)=[];
%% Plotting
% figure;
% imagesc(P(2:end-1,2:end-1,2));
% axes1 = gca;
% box(axes3,'on');
% set(axes3,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
% title('P')
% axes3.YDir='normal';
% colorbar
% save('Re1000Nodes200ErrorE5.mat','u','v','P','Error2','dt','dx','P0','Re')
uCentral=interp2(uXlocations,uYlocations,u(:,:,k),pXlocations,pYlocations);
vCentral=interp2(vXlocations,vYlocations,v(:,:,k),pXlocations,pYlocations);

figure;
imagesc(uCentral(2:end-1,2:end-1));
axes1 = gca;
box(axes1,'on');
set(axes1,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
title('u')
axes1.YDir='normal';
colorbar
yPlot=1:Nodes;xPlot=1:Nodes;
% Horizontal grid
for k = 1:length(yPlot)
line([xPlot(1) xPlot(end)], [yPlot(k) yPlot(k)],'Color',[0 0 0],'LineWidth',1)
end
% Vertical grid
for k = 1:length(xPlot)
line([xPlot(k) xPlot(k)], [yPlot(1) yPlot(end)],'Color',[0 0 0],'LineWidth',1)
end
xlim([0 Nodes-1])
ylim([0 Nodes-1])
% NodeInterp=
NodeLabelsNum=1:6:Nodes;
NodeLabels=num2cell(NodeLabelsNum);
for i=1:length(NodeLabels)
NodeLabels{i}=num2str(NodeLabels{i});
end
set(axes1,'FontName','Times New Roman','FontSize',25,'Layer','top','LineWidth',3,'XTick',NodeLabelsNum,...
    'XTickLabel',NodeLabels,'YTick',NodeLabelsNum,'YTickLabel',NodeLabels);
xlabel('x')
ylabel('y')
axes1.TickLength=[0 0];

figure;
imagesc(vCentral(2:end-1,2:end-1));
axes2 = gca;
box(axes2,'on');
set(axes2,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
colorbar
title('v')
axes2.YDir='normal';
yPlot=1:Nodes;xPlot=1:Nodes;
% Horizontal grid
for k = 1:length(yPlot)
line([xPlot(1) xPlot(end)], [yPlot(k) yPlot(k)],'Color',[0 0 0],'LineWidth',1)
end
% Vertical grid
for k = 1:length(xPlot)
line([xPlot(k) xPlot(k)], [yPlot(1) yPlot(end)],'Color',[0 0 0],'LineWidth',1)
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
