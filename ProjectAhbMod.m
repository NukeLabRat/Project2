%% CFD PROJECT 2
% Dan Gould

close all
clear
clc

Re=400;
dt=.008;
TimeSteps=1;
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
PoissonIn.dx=dx;
PoissonIn.dy=dy;
PoissonIn.xSize=pSize(2);
PoissonIn.ySize=pSize(1);

IsCenterX=logical(IsCenterX);
IsCenterY=logical(IsCenterY);
IsCenterP=logical(IsCenterP);
u=zeros(uSize(1),uSize(2),TimeSteps);
v=zeros(vSize(1),vSize(2),TimeSteps);
P=ones(pSize(1),pSize(2),TimeSteps);
ConstantMat=zeros(pSize(1),pSize(2));
Ustar(:,:,1)=zeros(uSize(1),uSize(2));
Vstar(:,:,1)=zeros(vSize(1),vSize(2));
p=P(:,:,1);
C=p;

% for i = 1:xEnd-1 %Assign velocity BC for initial Step
%     for j = 1:yEnd
%         if IsCenterX(j,i)==true %checks if node is central node
%         else %For Boundary Nodes
%             if j==1
%                 u(j,i)=-u(2,i)/2;
%             elseif j==uSize(1)
%                 u(j,i)=2-u(uSize(1)-1,i);
%             else
%                 u(j,i)=0;
%             end
%         end
%     end
% end
PoissonError=1;
StartingTime=tic;
TimeCheck=0;
Error2=1;
MainIterations=1;
while Error2>5E-7 || MainIterations<100
    
    PoissonIn.PoissonErrorMax=Error2/50;
    if Error2<.05
        PoissonIn.SOR=1.5;
    else
        PoissonIn.SOR=1;
    end

    v(:,1)=-v(:,2);
    v(:,end)=-v(:,end-1);
    v(1,:)=0;
    v(end,:)=0;
    
    u(:,1)=0;
    u(:,end)=0;
    u(end,:)=2-u(end-1,:);
    u(1,:)=-u(2,:);
    
    
    P(:,:,1)=P(:,:,1);
    
    unow=u(:,:);
    vnow=v(:,:);
    pnow=P(:,:);
    uCentralStart=interp2(uXlocations,uYlocations,u,pXlocations,pYlocations);
    vCentralStart=interp2(vXlocations,vYlocations,v,pXlocations,pYlocations);
    SpeedOld=(uCentralStart(IsCenterP).^2+vCentralStart(IsCenterP).^2).^.5;
    
    
    for j = 2:yEnd-1
        for i = 2:xEnd-2
%             if IsCenterX(j,i)==true %checks if node is central node
                Ustar(j,i)=u(j,i)+((dt/(Re*dx^2))*(u(j,i+1)-2*u(j,i)+u(j,i-1)))...
                    +(((dt/(Re*dy^2)))*(u(j+1,i)-2*u(j,i)+u(j-1,i)))...
                    -(dt/dx)*(((u(j,i)+u(j,i+1))/2)^2-((u(j,i-1)+u(j,i))/2)^2)...
                    -(dt/dy)*(((v(j,i)+v(j,i+1))/2)*((u(j,i)+u(j+1,i))/2)-(((v(j-1,i)+v(j,i+1))/2)*((u(j-1,i)+u(j,i))/2)));
%             end
        end
    end
    
    Ustar(end,:)=2-Ustar(end-1,:);
    Ustar(1,:)=-Ustar(2,:);
    Ustar(:,1)=0.0;
    Ustar(:,end)=0.0;
    
     for j = 2:yEnd-2
        for i = 2:xEnd-1
%             if IsCenterY(j,i)==true %checks if node is central node
                Vstar(j,i)=v(j,i)+((dt/(Re*dx^2))*(v(j,i+1)-2*v(j,i)+v(j,i-1)))...
                    +(((dt/(Re*dy^2)))*(v(j+1,i)-2*v(j,i)+v(j-1,i)))...
                    -(dt/dy)*(((v(j,i)+v(j+1,i))/2)^2-((v(j-1,i)+v(j,i))/2)^2)...
                    -(dt/dx)*((((u(j,i)+u(j+1,i))/2)*((v(j,i)+v(j,i+1))/2))-(((u(j,i-1)+u(j+1,i-1))/2)*((v(j,i-1)+v(j,i))/2)));
                
                %             else %For Boundary Nodes
                %                 Vstar(j,i)=0;
%             end
        end
    end
    Vstar(1,:)=0.0;
    Vstar(end,:)=0.0;
    Vstar(:,1)=-Vstar(:,2);
    Vstar(:,end)=-Vstar(:,end-1);
    
%     ConstantMat=padarray((diff(Ustar(2:end-1,:),1,2)/dx+diff(Vstar(:,2:end-1),1,1)/dy)./dt,[1 1]);
%     P0=P(:,:);
    PoissonError=1;
    while PoissonError > PoissonIn.PoissonErrorMax
        p(:,1) = p(:,2);
        p(:,end) = p(:,end-1);
        p(1,:) = p(2,:);
        p(end,:) = p(end-1,:);
        p1 = p;
        
        for i=2:xEnd-1
            
            for j=2:yEnd-1
                C(j,i)=((Ustar(j,i)-Ustar(j,i-1))/(dx*dt))+((Vstar(j,i)-Vstar(j-1,i))/(dy*dt));
                p2(j,i) = (1-PoissonIn.SOR)*p(j,i)+PoissonIn.SOR*(((dy.^2)*(p(j,i-1)+p(j,i+1))+(dx.^2)*(p(j-1,i)+p(j+1,i))-(dx.^2)*(dy.^2)*C(j,i))/(2*(dx.^2+dy.^2)));
                p(j,i) = p2(j,i);
            end
        end
        PoissonError = norm(p1-p,'fro');
    end
   
    %     [Pressure ~]=PoisonPressure6(ConstantMat,IsCenterP,P0,dx,dy,PoissonIn);
    
    
    %     norm(Pressure-p,'fro')
    Pressure=p;
    P(:,:,1)=p;
    for j=2:yEnd-1
        for i=2:xEnd-2
            u(j,i)=Ustar(j,i)-(dt/dx)*(p(j,i+1)-p(j,i));
        end
    end
    
    for j=2:yEnd-2
        for i=2:xEnd-1
            v(j,i)=Vstar(j,i)-(dt/dy)*(p(j+1,i)-p(j,i));
        end
    end
    %     for i = 1:xEnd-1
    %         for j = 1:yEnd
    %             if IsCenterX(j,i)==true %checks if node is central node
    %                 u(j,i,k+1)=Ustar(j,i)-(dt/dx)*(Pressure(j,i+1)-Pressure(j,i));
    %                 %             else %For Boundary Nodes
    %                 %                 if j==1
    %                 %                     u(j,i,k+1)=-u(2,i,k);
    %                 %                 elseif j==uSize(1)
    %                 %                     u(j,i,k+1)=2-u(uSize(1)-1,i,k);
    %                 %                 else
    %                 %                     u(j,i,k+1)=0;
    %                 %                 end
    %                 %
    %             end
    %
    %         end
    %     end
    %     for i = 1:xEnd
    %         for j = 1:yEnd-1
    %             if IsCenterY(j,i)==true %checks if node is central node
    %                 v(j,i,k+1)=Vstar(j,i)-(dt/dy)*(Pressure(j+1,i)-Pressure(j,i));
    %                 %             else %For Boundary Nodes
    %                 %                 if i==1
    %                 %                     v(j,i,k+1)=-v(j,2,k);
    %                 %                 elseif i==vSize(2)
    %                 %                     v(j,i,k+1)=-v(j,vSize(2)-1,k);
    %                 %                 else
    %                 %                     v(j,i,k+1)=0;
    %                 %                 end
    %             end
    %         end
    %     end
    %
    
%     v(:,1,2)=-v(:,2,2);
%     v(:,end,2)=-v(:,end-1,2);
%     v(1,:,2)=0;
%     v(end,:,2)=0;
%     
%     u(:,1,2)=0;
%     u(:,end,2)=0;
%     u(end,:,2)=2-u(end-1,:,2);
%     u(1,:,2)=-u(2,:,2);
    
    uCentralEnd=interp2(uXlocations,uYlocations,u(:,:),pXlocations,pYlocations);
    vCentralEnd=interp2(vXlocations,vYlocations,v(:,:),pXlocations,pYlocations);
%     
%   for j=1:yEnd-1
%     	for i=1:xEnd-1
%     		un(j,i)=(u(j,i)+u(j+1,i))/2;
%         end
%   end
%     for j=1:yEnd-1
%     	for i=1:xEnd-1
%     		vn(j,i)=((v(j,i)+v(j,i+1)))/2;
%         end
%     end
    
    
    SpeedNew=(uCentralEnd(IsCenterP).^2+vCentralEnd(IsCenterP).^2).^.5;
    Error2=norm(SpeedNew-SpeedOld,'fro');
    MainIterations=MainIterations+1;
    
    TimeCheck=TimeCheck+1;
    if TimeCheck==20
        clc
        Stab1
        Stab2
        CurrentTime=dt*MainIterations
        Error=Error2
        TimeCheck=0;
    end
end

ComputationTime=toc(StartingTime)
%% Plotting
figure;
imagesc(P(2:end-1,2:end-1));
axes1 = gca;
box(axes1,'on');
set(axes1,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
title('P')
axes1.YDir='normal';
colorbar

figure;
imagesc(u(2:end-1,:));
axes1 = gca;
box(axes1,'on');
set(axes1,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
title('u')
axes1.YDir='normal';
colorbar

figure;
imagesc(v(:,2:end-1));
axes1 = gca;
box(axes1,'on');
set(axes1,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
colorbar
title('v')
axes1.YDir='normal';

figure;
imagesc((uCentralEnd(2:end-1,2:end-1).^2+vCentralEnd(2:end-1,2:end-1).^2).^.5);
axes4 = gca;
box(axes4,'on');
set(axes4,'FontName','Times New Roman','FontSize',25,'LineWidth',3)
colorbar
title('speed')
axes4.YDir='normal';
