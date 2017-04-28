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


% u=zeros(yEnd,xEnd)+.25;
u=X;

uCentral=interp2(uXlocations,uYlocations,u,pXlocations,pYlocations);
vCentral=interp2(vXlocations,vYlocations,u,pXlocations,pYlocations);

P = zeros(yEnd-1,xEnd-1); %initialize matrix with zeros
Nodes = 1:numel(P);
%% Apply BC
Pinitial = 0;
BoundaryNodesP=MatEdges(P);
P(BoundaryNodesP) = Pinitial; %Assigns value of 1 to borders of matrix;

%% Center Nodes
CenterNodes = Nodes(~ismember(Nodes,BoundaryNodesP));

rng('default')
u=rand(yEnd,xEnd);
v=u;
Ustar=u;
Vstar=u;

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
            
            
            
            Ustar(j,i) = (-dxu2-uvdy(j,i)+1./Re.*(dudx2+dudy2)).*dt+u(j,i);
            Vstar(j,i) = (-dyv2-uvdx(j,i)+1./Re.*(dvdx2+dvdy2)).*dt+v(j,i);
            
        else %For Boundary Nodes
            Ustar(j,i) = 0;
            Vstar(j,i) =0;
            
        end
        
    end
end

uStarCentral=interp2(uXlocations,uYlocations,Ustar,pXlocations,pYlocations);
vStarCentral=interp2(vXlocations,vYlocations,Vstar,pXlocations,pYlocations);

for i = 1:xEnd
    for j = 1:yEnd
        if IsCenter(j,i)==true %checks if node is central node
            
            dxUstarCentral(j,i)=(uStarCentral(j,i+1)^2-uStarCentral(j,i)^2)/dx; %is in cell center
            dyVstarCentral(j,i)=(vStarCentral(j+1,i)^2-vStarCentral(j,i)^2)/dy;

            
         
            
        else %For Boundary Nodes
            dxUstarCentral(j,i)=0;
            dyVstarCentral(j,i)=0;
        end
        
    end
end





Ustar=3