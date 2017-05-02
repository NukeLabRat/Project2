function [ Pressure, Iterations] = PoisonPressureGPU( ConstantMat, IsCenterP, P0, xSize, ySize)
%PoisonPressure Pressure solving function
%   Itteratively solves for the pressure field durring each timestep. Gives
%   back the pressure field in a matrix at locations given in NodeX and
%   NodeY.

dx=.05;
dy=.05;
Iterations = 0;
Error2 = 1;
SOR=1.3; %1.7189 is optimal value.
Beta=dx/dy;
% [ySize, xSize] = size(IsCenterP);
% ConstantMat(isnan(ConstantMat))=0;
Pold=P0;
Pressure=Pold;
while Error2>1E-8
    for i = 1:xSize
        for j = 1:ySize
            if IsCenterP(j,i)==false %checks if node is central node
                if j==1
                    Pressure(j,i)=Pold(j+1,i);
%                     Pold(j,i)=Pold(j+1,i);
                end
                if j==ySize
                    Pressure(j,i)=Pold(j-1,i);
%                     Pold(j,i)=Pold(j-1,i);
                end
                if i==1
                    Pressure(j,i)=Pold(j,i+1);
%                     Pold(j,i)=Pold(j,i+1);
                end
                if i==xSize
                    Pressure(j,i)=Pold(j,i-1);
%                     Pold(j,i)=Pold(j,i-1);
                end
            end   
        end 
    end

    for i = (1:xSize)
        for j = (1:ySize)
            if IsCenterP(j,i)==true %checks if node is central node
%                  RHS=P(j,i-1)/dx^2-ConstantMat(j,i)+(P(j-1,i)+P(j+1,i))/dy^2+P(j,i+1)/dx^2;
%                  P(j,i)=RHS/(1/dx^2+2/dy^2+1/dx^2);
                 Pressure(j,i) = (1-SOR).*Pold(j,i)+SOR.*(Pold(j,i+1)+Pressure(j,i-1)+Beta.^2.*(Pold(j+1,i)+Pressure(j-1,i))+dx.^2.*ConstantMat(j,i))./(2.*(1+Beta.^2));
                
            end
        end
    end
    
    Error2 = norm(Pressure-Pold); %Calculate norm 2 error
    if Iterations ==50000
        Stop=1;
    end
    Pold=Pressure;
    Iterations = Iterations+1;
end
Done=1;
end

