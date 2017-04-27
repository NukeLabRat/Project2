%% CFD PROJECT 2

% Dan Gould

% close all
clear
clc

Re=1;
dt=.5;
%% Geometry - 
L = 4; %m, y-dir 
W = 6; %m, x-dir
C = 1.5; %m
A = 1.75; %m
B = 4.25; %m
D = 4.5; %m

%%
dx = 0.25; %m
dy = 0.25; %m
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
[uXlocations, uYlocations] = meshgrid(x,yPlusHalf);
[vXlocations, vYlocations] = meshgrid(xPlusHalf,y);
[uvXlocations, uvYlocations] = 

% u=zeros(yEnd,xEnd)+.25;
u=X;

uCentral=interp2(uXlocations,uYlocations,u,pXlocations,pYlocations);
vCentral=interp2(vXlocations,vYlocations,u,pXlocations,pYlocations);
Psi = zeros(yEnd,xEnd)+.5; %initialize matrix with zeros
Nodes = 1:numel(Psi);
%% Apply BC
BoundaryStream = 1;
BoundaryNodes=MatEdges(Psi);
[BoundaryY, BoundaryX] = ind2sub(size(Psi),BoundaryNodes);
Psi(BoundaryNodes) = BoundaryStream; %Assigns value of 1 to borders of matrix;

Slot1 = x>C & x<A; %determine indices of boundary nodes in slot 1
Slot2 = x>B & x<=D; %determine indices of boundary nodes in slot 2
Psi(1,Slot1)=interp1([C A],[BoundaryStream 0],x(Slot1));
Psi(1,Slot2)=interp1([B D],[0 BoundaryStream],x(Slot2));

Wall1 = x>=A & x<B; %determine indices of boundary nodes in between slots
Psi(1,Wall1) = 0;

%% Center Nodes
CenterNodes = Nodes(~ismember(Nodes,BoundaryNodes));
[CenterY, CenterX] = ind2sub(size(Psi),CenterNodes);

u=zeros(yEnd,xEnd);
v=u;
uStar=u;
vStar=u;

Index=zeros(yEnd,xEnd);
IsCenter=zeros(yEnd,xEnd);
    for i = 1:xEnd %This loop determines whether each node is central node or boundary node.
        for j = 1:yEnd
            Index(j,i) = sub2ind([yEnd xEnd],j,i);
            if sum(ismember(CenterNodes,Index(j,i)))==1
                IsCenter(j,i) = true;
            else
                IsCenter(j,i) = false;
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
        for j = 1:yEnd
            if IsCenter(j,i)==true %checks if node is central node
                                
                dxu2=(uCentral(j,i+1)^2-uCentral(j,i)^2)/dx;
                dyv2=(vCentral(j+1,i)^2-vCentral(j,i)^2)/dy;
                
                dudx2=(u(j,i+1)-2*u(j,i)+u(j,i-1))./dx^2;
                dudy2=(u(j-1,i)-2*u(j,i)+u(j+1,i))./dy^2;
                
                dvdx2=(v(j,i+1)-2*v(j,i)+v(j,i-1))./dx^2;
                dvdy2=(v(j-1,i)-2*v(j,i)+v(j+1,i))./dy^2;
                
                
                
                uStar(j,i) = (-dxu2-uvdy(j,i)+1./Re.*(dudx2+dudy2)).*dt+u(j,i);
                vStar(j,i) = (-dyv2-uvdx(j,i)+1./Re.*(dvdx2+dvdy2)).*dt+v(j,i);
                
                
                
                
            else
                Psi(j,i,k) = Psi(j,i,k-1); %applies BC condition of Boundary Nodes
            end
            
        end
    end
    
    
    
    
    
uStar=3