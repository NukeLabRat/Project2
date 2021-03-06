%% CFD PROJECT 2

% Dan Gould

 close all
clear
clc

Re=1;
dt=.0005;
TimeSteps=20;
%% Geometry -
L = 1; %m, y-dir
W = 1; %m, x-dir
%%
dx = 0.05; %m
dy = 0.05; %m
Beta = dx/dy;

x = 0:dx:W; %vector of node locations; x-dir
y = 0:dy:L; %vector of node locations; y-dir
xPlusHalf = x+.5*dx;
xMinusHalf = x-.5*dx;
yPlusHalf = y+.5*dy;
yMinusHalf = y-.5*dy;


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

NodesP = 1:numel(P);
NodesX = 1:numel(uXlocations);
NodesY = 1:numel(uXlocations);
%% Apply BC
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



u=zeros(uSize(1),uSize(2),TimeSteps);
v=zeros(vSize(1),vSize(2),TimeSteps);
P=zeros(pSize(1),pSize(2),TimeSteps);
P0=P(:,:,1);
Ustar(:,:,1)=ones(uSize(1),uSize(2));
Vstar(:,:,1)=ones(vSize(1),vSize(2));
dxUstar(:,:,1)=ones(uSize(1),uSize(2));
dyVstar(:,:,1)=ones(vSize(1),vSize(2));
duStarCentral(:,:,1)=ones(pSize(1),pSize(2));
dvStarCentral(:,:,1)=ones(pSize(1),pSize(2));
SOR=1.3;

    for i = 1:xEnd-1
        for j = 1:yEnd
            if IsCenterX(j,i)==true %checks if node is central node
                u(j,i,1) = 0;
            else %For Boundary Nodes
                if j==1
                    u(j,i,1)=0;
                elseif j==uSize(1)
                    u(j,i,1)=2;
                    u(j,i,2)=2;
                else
                    u(j,i)=0;
                end
                
            end
            
        end
    end

for k=2:TimeSteps
    
    uCentral=interp2(uXlocations,uYlocations,u(:,:,k-1),pXlocations,pYlocations);
    vCentral=interp2(vXlocations,vYlocations,v(:,:,k-1),pXlocations,pYlocations);
    uvdy=(interp2(uXlocations,uYlocations,u(:,:,k-1),uXlocations,uYlocations+.5*dy).*interp2(vXlocations,vYlocations,v(:,:,k-1),uXlocations,uYlocations+.5*dy)...
        -interp2(uXlocations,uYlocations,u(:,:,k-1),uXlocations,uYlocations-.5*dy).*interp2(vXlocations,vYlocations,v(:,:,k-1),uXlocations,uYlocations-.5*dy))/dy;
    uvdx=(interp2(uXlocations,uYlocations,u(:,:,k-1),vXlocations+.5*dx,vYlocations).*interp2(vXlocations,vYlocations,v(:,:,k-1),vXlocations+.5*dx,vYlocations)...
        -interp2(uXlocations,uYlocations,u(:,:,k-1),vXlocations-.5*dx,vYlocations).*interp2(vXlocations,vYlocations,v(:,:,k-1),vXlocations-.5*dx,vYlocations))/dx;
    
    for i = 1:xEnd-1
        for j = 1:yEnd
            if IsCenterX(j,i)==true %checks if node is central node

                dxu2=(uCentral(j,i+1)^2-uCentral(j,i)^2)/dx;
                dudx2=(u(j,i+1,k-1)-2*u(j,i,k-1)+u(j,i-1,k-1))./dx^2;
                dudy2=(u(j-1,i,k-1)-2*u(j,i,k-1)+u(j+1,i,k-1))./dy^2;
                Ustar(j,i) = (-dxu2-uvdy(j,i)+1./Re.*(dudx2+dudy2)).*dt+u(j,i,k-1);
            else %For Boundary Nodes
                Ustar(j,i)=0;
            end
            
        end
    end
    for i = 1:xEnd
        for j = 1:yEnd-1
            if IsCenterY(j,i)==true %checks if node is central node
                dyv2=(vCentral(j+1,i)^2-vCentral(j,i)^2)/dy;
                dvdx2=(v(j,i+1,k-1)-2*v(j,i,k-1)+v(j,i-1,k-1))./dx^2;
                dvdy2=(v(j-1,i,k-1)-2*v(j,i,k-1)+v(j+1,i,k-1))./dy^2;
                Vstar(j,i) = (-dyv2-uvdx(j,i)+1./Re.*(dvdx2+dvdy2)).*dt+v(j,i,k-1);
            else %For Boundary Nodes
                Vstar(j,i) =0;
            end
            
        end
    end
    
    for i = 1:xEnd-1
        for j = 1:yEnd
            if IsCenterX(j,i)==true %checks if node is central node
                dxUstar(j,i)=(Ustar(j,i+1)^2-Ustar(j,i)^2)/dx; %is in cell center
            else %For Boundary Nodes
                dxUstar(j,i)=0;
            end
        end
    end
    
    for i = 1:xEnd
        for j = 1:yEnd-1
            if IsCenterY(j,i)==true %checks if node is central node
                dyVstar(j,i)=(Vstar(j+1,i)^2-Vstar(j,i)^2)/dy;
            else %For Boundary Nodes
                dyVstar(j,i)=0;
            end
        end
    end
    

    
    duStarCentral=interp2(uXlocations,uYlocations,dxUstar,pXlocations,pYlocations);
    dvStarCentral=interp2(vXlocations,vYlocations,dyVstar,pXlocations,pYlocations);
    ConstantMat=(duStarCentral+dvStarCentral)./dt;
    P(:,:,k)=PoisonPressure(ConstantMat,IsCenterP,P0,dx,dy,SOR);
    P0=P(:,:,k);
    Px=diff(P(:,:,k))./dx;
    Py=diff(P(:,:,k)')./dy;
    PinterpU=interp2(pXlocations,pYlocations,P(:,:,k),uXlocations,uYlocations);
    PinterpV=interp2(pXlocations,pYlocations,P(:,:,k),vXlocations,vYlocations);
    for i = 1:xEnd-1
        for j = 1:yEnd
            if IsCenterX(j,i)==true %checks if node is central node
                u(j,i,k) = Ustar(j,i)+dt/dx*Px(j,i);
            else %For Boundary Nodes
                if j==1
                    u(j,i,k)=-u(2,i,k-1);
                elseif j==uSize(1)
                    u(j,i,k)=2-u(uSize(1)-1,i,k-1);
                else
                    u(j,i,k)=0;
                end
                
            end
            
        end
    end
    for i = 1:xEnd
        for j = 1:yEnd-1
            if IsCenterY(j,i)==true %checks if node is central node
                v(j,i,k) = Vstar(j,i)+dt/dy*Py(j,i);
            else %For Boundary Nodes
                if i==1
                    v(j,i,k)=-v(j,2,k-1);
                elseif i==vSize(2)
                    v(j,i,k)=-v(j,vSize(2)-1,k-1);
                else
                    v(j,i,k)=0;
                end
                

            end
            
        end
    end
    u(:,:,k+1)=u(:,:,k);
    v(:,:,k+1)=v(:,:,k);
    P(:,:,k+1)=P0;
end
%% Plotting
figure;
imagesc(P(2:end-1,2:end-1,k));
axes1 = gca;
box(axes1,'on');
set(axes1,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
title('P')
axes1.YDir='normal';
colorbar

figure;
imagesc(u(2:end-1,:,k));
axes1 = gca;
box(axes1,'on');
set(axes1,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
title('u')
axes1.YDir='normal';
colorbar

figure;
imagesc(v(:,2:end-1,k));
axes1 = gca;
box(axes1,'on');
set(axes1,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
colorbar
title('v')
axes1.YDir='normal';
