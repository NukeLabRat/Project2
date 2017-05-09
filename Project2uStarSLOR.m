%% CFD PROJECT 2
% Dan Gould

close all
clear
clc

Re=100;
dt=.0005;
TimeSteps=2;
Nodes=50;
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
PoissonIn.dx=dx;
PoissonIn.dy=dy;
PoissonIn.xSize=pSize(2);
PoissonIn.ySize=pSize(1);
IsCenterP=logical(IsCenterP);

PoissonIn.IsCenterP=logical(IsCenterP);
u=zeros(uSize(1),uSize(2),TimeSteps);
v=zeros(vSize(1),vSize(2),TimeSteps);
P=zeros(pSize(1),pSize(2),TimeSteps);
ConstantMat=zeros(pSize(1),pSize(2));
% ConstantMat2=zeros(pSize(1),pSize(2));
Ustar(:,:,1)=zeros(uSize(1),uSize(2));
Vstar(:,:,1)=zeros(vSize(1),vSize(2));
PressureIterations=zeros(1,10000);


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

StartingTime=tic;
TimeCheck=0;
Error2=1;
MainIterations=1;
while Error2>5E-7 || MainIterations<100
PoissonIn.PoissonErrorMax=Error2/50;
    
    
    if Error2<.05
        PoissonIn.SOR=1.2;
    else
        PoissonIn.SOR=1;
    end
    

    u(:,:,1)=u(:,:,2);
    v(:,:,1)=v(:,:,2);
    P(:,:,1)=P(:,:,2);
    k=1;
    unow=u(:,:,1);
    vnow=v(:,:,1);
    pnow=P(:,:,1);
%     uCentral=interp2(uXlocations,uYlocations,u(:,:,k),pXlocations,pYlocations);
%     vCentral=interp2(vXlocations,vYlocations,v(:,:,k),pXlocations,pYlocations);
% %     uvdy=(interp2(uXlocations,uYlocations,u(:,:,k),uXlocations,uYlocations+.5*dy).*interp2(vXlocations,vYlocations,v(:,:,k),uXlocations,uYlocations+.5*dy)...
%         -interp2(uXlocations,uYlocations,u(:,:,k),uXlocations,uYlocations-.5*dy).*interp2(vXlocations,vYlocations,v(:,:,k),uXlocations,uYlocations-.5*dy))/dy;
%     uvdx=(interp2(uXlocations,uYlocations,u(:,:,k),vXlocations+.5*dx,vYlocations).*interp2(vXlocations,vYlocations,v(:,:,k),vXlocations+.5*dx,vYlocations)...
%         -interp2(uXlocations,uYlocations,u(:,:,k),vXlocations-.5*dx,vYlocations).*interp2(vXlocations,vYlocations,v(:,:,k),vXlocations-.5*dx,vYlocations))/dx;
%     dxuSquared=(uCentral.^2-uCentral.^2)/dx;
%     dyvSquared=(vCentral.^2-vCentral.^2)/dy;
    
    for i = 1:xEnd-1
        for j = 1:yEnd
            if IsCenterX(j,i)==true %checks if node is central node
%                 du2dx2=(u(j,i+1,k)-2*u(j,i,k)+u(j,i-1,k))./dx^2;
%                 du2dy2=(u(j+1,i,k)-2*u(j,i,k)+u(j-1,i,k))./dy^2;
%                 Ustar(j,i) = (-dxuSquared(j,i)-uvdy(j,i)+1./Re.*(du2dx2+du2dy2)).*dt+u(j,i,k);
                
                     	Ustar(j,i)=u(j,i,k)+((dt/(Re*dx^2))*(u(j,i+1,k)-2*u(j,i,k)+u(j,i-1,k)))...
                +(((dt/(Re*dy^2)))*(u(j+1,i,k)-2*u(j,i,k)+u(j-1,i,k)))...
                -(dt/dx)*(((u(j,i,k)+u(j,i+1,k))/2)^2-((u(j,i-1,k)+u(j,i,k))/2)^2)...
                -(dt/dy)*(((v(j,i,k)+v(j,i+1,k))/2)*((u(j,i,k)+u(j+1,i,k))/2)-(((v(j-1,i,k)+v(j,i+1,k))/2)*((u(j-1,i,k)+u(j,i,k))/2)));
                
                
            else %For Boundary Nodes
%                 if j==uSize(1)
%                     Ustar(j,i)=2-Ustar(yEnd-1,i);
% %                     Ustar2(j,i)=2-Ustar2(yEnd-1,i,k);
%                 else
                    Ustar(j,i)=0;
% %                     Ustar2(j,i)=0;
%                 end
                
                %                 Ustar(j,i)=0;
            end
            
        end
    end
    
    Ustar(1,:)=-Ustar(2,:);
    
    Ustar(uSize(1),:)=2-Ustar(yEnd-1,:);
%       Ustar2(1,:)=-Ustar2(2,:);
%     
%     Ustar2(uSize(1),:)=2-Ustar2(uSize(1)-1,:);
    for i = 1:xEnd
        for j = 1:yEnd-1
            if IsCenterY(j,i)==true %checks if node is central node
%                 dv2dx2=(v(j,i+1,k)-2*v(j,i,k)+v(j,i-1,k))./dx^2;
%                 dv2dy2=(v(j+1,i,k)-2*v(j,i,k)+v(j-1,i,k))./dy^2;
%                 Vstar(j,i) = (-dyvSquared(j,i)-uvdx(j,i)+1./Re.*(dv2dx2+dv2dy2)).*dt+v(j,i,k);
                
                        	Vstar(j,i)=v(j,i,k)+((dt/(Re*dx^2))*(v(j,i+1,k)-2*v(j,i,k)+v(j,i-1,k)))...
                +(((dt/(Re*dy^2)))*(v(j+1,i,k)-2*v(j,i,k)+v(j-1,i,k)))...
                -(dt/dy)*(((v(j,i,k)+v(j+1,i,k))/2)^2-((v(j-1,i,k)+v(j,i,k))/2)^2)...
                -(dt/dx)*((((u(j,i,k)+u(j+1,i,k))/2)*((v(j,i,k)+v(j,i+1,k))/2))-(((u(j,i-1,k)+u(j+1,i-1,k))/2)*((v(j,i-1,k)+v(j,i,k))/2)));

            else %For Boundary Nodes
%                 Vstar(j,i)=0;
            end
        end
    end
        Vstar(1,:)=0.0;
    Vstar(end,:)=0.0;
    Vstar(:,1,k)=-Vstar(:,2);
    Vstar(:,end)=-Vstar(:,end-1);
    
%     Vstar(:,1)=-Vstar(:,end);
%     Vstar(:,1)=-Vstar(:,2);
%     
%     Vstar(:,vSize(2))=-Vstar(:,vSize(2)-1);
%     
%     VstarNorm=norm(Vstar-Vstar2,inf)
    PoissonIn.ConstantMat=padarray((diff(Ustar(2:end-1,:),1,2)/dx+diff(Vstar(:,2:end-1),1,1)/dy)./dt,[1 1]);
    PoissonIn.P0=P(:,:,k);
%     for i = 2:Nodes
%         for j = 2:Nodes
% 
%             ConstantMat2(j,i)=-((Ustar(j,i)-Ustar(j,i-1))/(dx*dt))+((Vstar(j,i)-Vstar(j-1,i))/(dy*dt));
%         end
%     end
    PoissonIn.Iterations=MainIterations;

    [Pressure, PoissonIterations]=PoisonPressureSLOR4(PoissonIn);
%     PressureIterations(MainIterations)=PoissonIterations;
    P(:,:,k+1)=Pressure;
    PinterpU=interp2(pXlocations,pYlocations,P(:,:,k+1),uXlocations,uYlocations);
    PinterpV=interp2(pXlocations,pYlocations,P(:,:,k+1),vXlocations,vYlocations);
    for i = 1:xEnd-1
        for j = 1:yEnd
            if IsCenterX(j,i)==true %checks if node is central node
%                 u(j,i,k+1) = Ustar(j,i)+dt/dx*(PinterpU(j,i+1)-PinterpU(j,i));
                u(j,i,k+1)=Ustar(j,i)-(dt/dx)*(Pressure(j,i+1)-Pressure(j,i));
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
%                 v(j,i,k+1) = Vstar(j,i)+dt/dy*(PinterpV(j+1,i)-PinterpV(j,i));
                v(j,i,k+1)=Vstar(j,i)-(dt/dy)*(Pressure(j+1,i)-Pressure(j,i));
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
    
    uCentralEnd=interp2(uXlocations,uYlocations,u(:,:,2),pXlocations,pYlocations);
    vCentralEnd=interp2(vXlocations,vYlocations,v(:,:,2),pXlocations,pYlocations);
    pnow=P(:,:,2);
    Error2=norm(u(:,:,2)-u(:,:,1),'fro');
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
uCentral=interp2(uXlocations,uYlocations,u(:,:,k),pXlocations,pYlocations);
vCentral=interp2(vXlocations,vYlocations,v(:,:,k),pXlocations,pYlocations);
ComputationTime=toc(StartingTime)
PressureIterations(PressureIterations==0)=[];
%% Plotting
figure;
imagesc(P(2:end-1,2:end-1,2));
axes1 = gca;
box(axes1,'on');
set(axes1,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
title('P')
axes1.YDir='normal';
colorbar

figure;
imagesc(u(2:end-1,:,2));
axes1 = gca;
box(axes1,'on');
set(axes1,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
title('u')
axes1.YDir='normal';
colorbar

figure;
imagesc(v(:,2:end-1,2));
axes1 = gca;
box(axes1,'on');
set(axes1,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
colorbar
title('v')
axes1.YDir='normal';
