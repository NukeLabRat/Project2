%% CFD PROJECT 2
% Dan Gould

close all
clear
clc

Re=200;
dt=.006;
TimeSteps=2;
Nodes=30;
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

IsCenterP=logical(IsCenterP);
u=zeros(uSize(1),uSize(2),TimeSteps);
v=zeros(vSize(1),vSize(2),TimeSteps);
P=zeros(pSize(1),pSize(2),TimeSteps);
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


Error2=1;
MainIterations=1;
while Error2>1E-4 || MainIterations<1000
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
                Ustar(j,i)=-u(j,i,k);
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
                Vstar(j,i) =-v(j,i,1);
%                 Vstar(j,i) =0;
            end
        end
    end
    
    for i = 1:xEnd-1
        for j = 1:yEnd
            if IsCenterX(j,i)==true %checks if node is central node
                dxUstar(j,i)=(Ustar(j,i+1)-Ustar(j,i))/dx; %is in cell center
            else %For Boundary Nodes
                dxUstar(j,i)=1;
            end
        end
    end
    
    for i = 1:xEnd
        for j = 1:yEnd-1
            if IsCenterY(j,i)==true %checks if node is central node
                dyVstar(j,i)=(Vstar(j+1,i)-Vstar(j,i))/dy;
            else %For Boundary Nodes
                dyVstar(j,i)=1;
            end
        end
    end
    ConstantMat=padarray((diff(Ustar(2:end-1,:),1,2)/dx+diff(Vstar(:,2:end-1),1,1)/dy)./dt,[1 1]);
    duStarCentral=interp2(uXlocations,uYlocations,dxUstar,pXlocations,pYlocations);
    dvStarCentral=interp2(vXlocations,vYlocations,dyVstar,pXlocations,pYlocations);
%     ConstantMat=(duStarCentral+dvStarCentral)./dt;

    [Pressure ~]=PoisonPressure3(ConstantMat,IsCenterP,P0,dx,dy);
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
    Error2=norm(P(:,:,2)-P(:,:,1));
    MainIterations=MainIterations+1;
%     P(:,:,k+1)=P(:,:,k);
%     CurrentTime=dt*MainIterations
end
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
