function [ Pressure, Iterations] = PoisonPressure6( ConstantMat, IsCenterP, P0, dx, dy, Data)
%PoisonPressure Pressure solving function
%   Itteratively solves for the pressure field durring each timestep. Gives
%   back the pressure field in a matrix at locations given in NodeX and
%   NodeY.

Iterations = 0;
Error2 = 1;
SOR=1; %1.7189 is optimal value.
Beta=dx/dy;
[ySize, xSize] = size(IsCenterP);
Pold=P0;
Pressure=Pold;
% Pressure2=Pold;
% Pold2=P0;
while Error2>Data.PoissonErrorMax
%     for i = 1:xSize
%         for j = 1:ySize
%             if IsCenterP(j,i)==false %checks if node is central node
%                 if j==1
%                     Pressure(j,i)=Pold(j+1,i);
%                     Pold(j,i)=Pold(j+1,i);
%                 end
%                 if j==ySize
%                     Pressure(j,i)=Pold(j-1,i);
%                     Pold(j,i)=Pold(j-1,i);
%                 end
%                 if i==1
%                     Pressure(j,i)=Pold(j,i+1);
%                     Pold(j,i)=Pold(j,i+1);
%                 end
%                 if i==xSize
%                     Pressure(j,i)=Pold(j,i-1);
%                     Pold(j,i)=Pold(j,i-1);
%                 end
%             end   
%         end 
%     end
        Pold(:,1) = Pold(:,2);
        Pold(:,end) = Pold(:,end-1);
        Pold(1,:) = Pold(2,:);
        Pold(end,:) = Pold(end-1,:);
        
%         Pold2(:,1) = Pold2(:,2);
%         Pold2(:,end) = Pold2(:,end-1);
%         Pold2(1,:) = Pold2(2,:);
%         Pold2(end,:) = Pold2(end-1,:);
%         P1=Pold;
%         P2=Pold2;

    for i = (1:xSize)
        for j = (1:ySize)
            if IsCenterP(j,i)==true %checks if node is central node
%                  RHS=P(j,i-1)/dx^2-ConstantMat(j,i)+(P(j-1,i)+P(j+1,i))/dy^2+P(j,i+1)/dx^2;
%                  P(j,i)=RHS/(1/dx^2+2/dy^2+1/dx^2);
                 Pressure(j,i) = (1-SOR).*Pold(j,i)+SOR.*(Pold(j,i+1)+Pold(j,i-1)+Beta^2.*(Pold(j+1,i)+Pold(j-1,i))-dx^2.*ConstantMat(j,i))./(2*(1+Beta^2));
%   Pressure2(j,i) = (1-SOR)*Pold2(j,i)+SOR*(((dy^2)*(Pold2(j,i-1)+Pold2(j,i+1))+(dx^2)*(Pold2(j-1,i)+Pold2(j+1,i))-(dx^2)*(dy^2)*ConstantMat(j,i))/(2*(dx^2+dy^2)));
  
            end
        end
    end
%     PressureDiff=norm(Pressure(2:end-1,2:end-1)-Pressure2(2:end-1,2:end-1),inf);
    Error2 = norm(Pressure(2:end-1,2:end-1)-Pold(2:end-1,2:end-1),'fro'); %Calculate norm 2 error
    if Iterations ==80000
        Stop=1;
    end
    if Iterations>100 && Error2>1E10
        Stop=1;
    end
    Pold=Pressure;
%     Pold2=Pressure2;
    Iterations = Iterations+1;
end
Done=1;
end

