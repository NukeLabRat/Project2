%% CFD PROJECT 2

% Dan Gould

% close all
clear
clc

Re=1;
dt=.5;
%% Geometry -
L = 1; %m, y-dir
W = 1; %m, x-dir

%%
dx = 0.1; %m
dy = 0.1; %m
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

rng('default')
u=rand(uSize(1),uSize(2));
v=u';
% Ustar=u;
% Vstar=u;


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


TimeSteps=10;
u=zeros(uSize(1),uSize(2),TimeSteps);
v=zeros(vSize(1),vSize(2),TimeSteps);
P=zeros(pSize(1),pSize(2),TimeSteps);
 Ustar(:,:,1)=ones(uSize(1),uSize(2));
 Vstar(:,:,1)=ones(vSize(1),vSize(2));
  dxUstarCentral(:,:,1)=ones(pSize(1),pSize(2));
 dyVstarCentral(:,:,1)=ones(pSize(1),pSize(2));

for k=2:TimeSteps
    
    u(1,:,k)=-u(2,:,k-1);
    u(uSize(1),:,k)=2-u(uSize(1)-1,:,k-1);
    v(:,2,k)=-v(:,2,k-1);
    v(:,vSize(2),k)=-v(:,vSize(2)-1,k-1);
    
    
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
                dudx2=(u(j,i+1,k)-2*u(j,i,k)+u(j,i-1,k))./dx^2;
                dudy2=(u(j-1,i,k)-2*u(j,i)+u(j+1,i,k))./dy^2;
                Ustar(j,i) = (-dxu2-uvdy(j,i)+1./Re.*(dudx2+dudy2)).*dt+u(j,i,k);
            else %For Boundary Nodes
                Ustar(j,i)=1;
            end
            
        end
    end
    for i = 1:xEnd
        for j = 1:yEnd-1
            if IsCenterY(j,i)==true %checks if node is central node
                dyv2=(vCentral(j+1,i)^2-vCentral(j,i)^2)/dy;
                dvdx2=(v(j,i+1,k)-2*v(j,i)+v(j,i-1,k))./dx^2;
                dvdy2=(v(j-1,i,k)-2*v(j,i)+v(j+1,i,k))./dy^2;
                Vstar(j,i) = (-dyv2-uvdx(j,i)+1./Re.*(dvdx2+dvdy2)).*dt+v(j,i,k);
            else %For Boundary Nodes
                
                Vstar(j,i) =1;
                
            end
            
        end
    end
    
    uStarCentral=interp2(uXlocations,uYlocations,Ustar,pXlocations,pYlocations);
    vStarCentral=interp2(vXlocations,vYlocations,Vstar,pXlocations,pYlocations);
    
    for i = 1:xEnd
        for j = 1:yEnd
            if IsCenterP(j,i)==true %checks if node is central node
                dxUstarCentral(j,i)=(uStarCentral(j,i+1)^2-uStarCentral(j,i)^2)/dx; %is in cell center
                dyVstarCentral(j,i)=(vStarCentral(j+1,i)^2-vStarCentral(j,i)^2)/dy;
                
            else %For Boundary Nodes
                dxUstarCentral(j,i)=0;
                dyVstarCentral(j,i)=0;
            end
        end
    end
    
    ConstantMat=(dxUstarCentral+dyVstarCentral)./dt;
    P(:,:,k)=PoisonPressure(ConstantMat,IsCenterP,P0);
    
end




Ustar=3