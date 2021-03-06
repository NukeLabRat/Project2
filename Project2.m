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
dx = 0.05; %m
dy = 0.05; %m
Beta = dx/dy;

x = 0:dx:W; %vector of node locations; x-dir
y = 0:dy:L; %vector of node locations; y-dir
xPlusHalf = x+.5*dx;
xMinusHalf = x-.5*dx;
yPlusHalf = y+.5*dy;
yMinusHalf = y-.5*dy;


xEnd = length(x); %length of location vectors
yEnd = length(y);

[X, Y] = meshgrid(x,y);
[pXlocations, pYlocations] = meshgrid(xPlusHalf(1:end-1),yPlusHalf(1:end-1));
[uXlocations, uYlocations] = meshgrid(x,yPlusHalf(1:end-1));
[vXlocations, vYlocations] = meshgrid(xPlusHalf(1:end-1),y);
uSize=size(uXlocations);
vSize=size(vXlocations);

P = zeros(yEnd-1,xEnd-1); %initialize matrix with zeros
P0=P;
NodesP = 1:numel(P);
NodesX = 1:numel(uXlocations);
NodesY = 1:numel(uXlocations);
%% Apply BC
Pinitial = 0;
BoundaryNodesP=MatEdges(P);
P(BoundaryNodesP) = Pinitial; %Assigns value of 1 to borders of matrix;
BoundaryNodesX=MatEdges(uXlocations);
BoundaryNodesY=MatEdges(vYlocations);

%% Center Nodes
CenterNodesP = NodesP(~ismember(NodesP,BoundaryNodesP));
CenterNodesX = NodesX(~ismember(NodesX,BoundaryNodesX));
CenterNodesY = NodesY(~ismember(NodesY,BoundaryNodesY));

rng('default')
u=rand(yEnd-1,xEnd);
v=u';
Ustar=u;
Vstar=u;
for TimeStep=1:10
    
    IndexX=zeros(yEnd-1,xEnd);
    IsCenterX=zeros(yEnd-1,xEnd);
    for i = 1:xEnd %This loop determines whether each node is central node or boundary node.
        for j = 1:yEnd-1
            IndexX(j,i) = sub2ind([yEnd-1 xEnd],j,i);
            if sum(ismember(CenterNodesX,IndexX(j,i)))==1
                IsCenterX(j,i) = true;
            else
                IsCenterX(j,i) = false;
            end
        end
    end
    IndexY=zeros(yEnd,xEnd-1);
    IsCenterY=zeros(yEnd,xEnd-1);
    for i = 1:xEnd-1 %This loop determines whether each node is central node or boundary node.
        for j = 1:yEnd
            IndexY(j,i) = sub2ind([yEnd xEnd-1],j,i);
            if sum(ismember(CenterNodesY,IndexY(j,i)))==1
                IsCenterY(j,i) = true;
            else
                IsCenterY(j,i) = false;
            end
        end
    end
    
    IndexP=zeros(yEnd-1,xEnd-1);
    IsCenterP=zeros(yEnd-1,xEnd-1);
    for i = 1:xEnd-1 %This loop determines whether each node is central node or boundary node.
        for j = 1:yEnd-1
            IndexP(j,i) = sub2ind([yEnd-1 xEnd-1],j,i);
            if sum(ismember(CenterNodesP,IndexP(j,i)))==1
                IsCenterP(j,i) = true;
            else
                IsCenterP(j,i) = false;
            end
        end
    end
    
    uCentral=interp2(uXlocations,uYlocations,u,pXlocations,pYlocations);
    vCentral=interp2(vXlocations,vYlocations,v,pXlocations,pYlocations);
    uvdy=(interp2(uXlocations,uYlocations,u,uXlocations,uYlocations+.5*dy).*interp2(vXlocations,vYlocations,v,uXlocations,uYlocations+.5*dy)...
        -interp2(uXlocations,uYlocations,u,uXlocations,uYlocations-.5*dy).*interp2(vXlocations,vYlocations,v,uXlocations,uYlocations-.5*dy))/dy;
    uvdx=(interp2(uXlocations,uYlocations,u,vXlocations+.5*dx,vYlocations).*interp2(vXlocations,vYlocations,v,vXlocations+.5*dx,vYlocations)...
        -interp2(uXlocations,uYlocations,u,vXlocations-.5*dx,vYlocations).*interp2(vXlocations,vYlocations,v,vXlocations-.5*dx,vYlocations))/dx;
    
    for i = 1:xEnd
        for j = 1:yEnd-1
            if IsCenterX(j,i)==true %checks if node is central node
                dxu2=(uCentral(j,i+1)^2-uCentral(j,i)^2)/dx;
                dudx2=(u(j,i+1)-2*u(j,i)+u(j,i-1))./dx^2;
                dudy2=(u(j-1,i)-2*u(j,i)+u(j+1,i))./dy^2;
                Ustar(j,i) = (-dxu2-uvdy(j,i)+1./Re.*(dudx2+dudy2)).*dt+u(j,i);
            else %For Boundary Nodes
                if i~=1 && i~=xEnd
                    if j==1
                     
                        dxu2=(uCentral(j,i+1)^2-uCentral(j,i)^2)/dx;
                        
                        dudx2=(u(j,i+1)-2*u(j,i)+u(j,i-1))./dx^2;
                        dudy2=((-u(j,i))-2*u(j,i)+(2-u(j,i)))./dy^2;
                        Ustar(j,i) = (-dxu2-uvdy(j,i)+1./Re.*(dudx2+dudy2)).*dt+u(j,i);
                    elseif j==yEnd-1
                        dxu2=(uCentral(j,i+1)^2-uCentral(j,i)^2)/dx;
                        dudx2=(u(j,i+1)-2*u(j,i)+u(j,i-1))./dx^2;
                        dudy2=(u(j-1,i)-2*u(j,i)+(2-u(j,i)))./dy^2;
                        Ustar(j,i) = (-dxu2-uvdy(j,i)+1./Re.*(dudx2+dudy2)).*dt+u(j,i);
                    else
                        Ustar(j,i)=0;
                    end
                else
                    Ustar(j,i)=0;
                end
                
            end
            
        end
    end
    for i = 1:xEnd-1
        for j = 1:yEnd
            if IsCenterY(j,i)==true %checks if node is central node
                dyv2=(vCentral(j+1,i)^2-vCentral(j,i)^2)/dy;
                dvdx2=(v(j,i+1)-2*v(j,i)+v(j,i-1))./dx^2;
                dvdy2=(v(j-1,i)-2*v(j,i)+v(j+1,i))./dy^2;
                Vstar(j,i) = (-dyv2-uvdx(j,i)+1./Re.*(dvdx2+dvdy2)).*dt+v(j,i);
            else %For Boundary Nodes
                if j~=1 && j~=yEnd
                    if i==1
                        dyv2=(vCentral(j+1,i)^2-vCentral(j,i)^2)/dy;
                        dvdx2=(v(j,i+1)-2*v(j,i)+(-v(j,i)))./dx^2;
                        dvdy2=(v(j-1,i)-2*v(j,i)+v(j+1,i))./dy^2;
                        Vstar(j,i) = (-dyv2-uvdx(j,i)+1./Re.*(dvdx2+dvdy2)).*dt+v(j,i);
                    elseif i==xEnd-1
                        dyv2=(vCentral(j+1,i)^2-vCentral(j,i)^2)/dy;
                        dvdx2=((-v(j,i))-2*v(j,i)+v(j,i-1))./dx^2;
                        dvdy2=(v(j-1,i)-2*v(j,i)+v(j+1,i))./dy^2;
                        Vstar(j,i) = (-dyv2-uvdx(j,i)+1./Re.*(dvdx2+dvdy2)).*dt+v(j,i);
                    else
                        Vstar(j,i) =0;
                    end
                else
                    Vstar(j,i) =0;
                end
            end
            
        end
    end
    
    uStarCentral=interp2(uXlocations,uYlocations,Ustar,pXlocations,pYlocations);
    vStarCentral=interp2(vXlocations,vYlocations,Vstar,pXlocations,pYlocations);
    
    for i = 1:xEnd-1
        for j = 1:yEnd-1
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
    P=PoisonPressure(ConstantMat,IsCenterP,P0);
end




Ustar=3