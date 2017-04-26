%% CFD PROJECT 2

% Dan Gould

% close all
clear
clc


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
[X, Y] = meshgrid(x,y);
xEnd = length(x); %length of location vectors
yEnd = length(y); 

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



uStar=