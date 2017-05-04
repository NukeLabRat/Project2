function [ Pressure, Iterations] = PoisonPressure5( ConstantMat, IsCenterP, P0, dx, dy,Data)
%PoisonPressure Pressure solving function
%   Itteratively solves for the pressure field durring each timestep. Gives
%   back the pressure field in a matrix at locations given in NodeX and
%   NodeY.

Iterations = 0;
Error2 = 1;
SOR=Data.SOR; 
Beta=dx/dy;
[ySize, xSize] = size(IsCenterP);
Pold=P0;
Pressure=Pold;
BetaSquared=Beta^2;
Divisor=2.*BetaSquared+2;
C=Data.ConstantMat;
while Error2>1E-8
    
%     Pold(Data.TopWallP)=Pold(Data.TopWallPmirror);
%     Pold(Data.BottomWallP)=Pold(Data.BottomWallPmirror);
%     Pold(Data.LeftWallP)=Pold(Data.LeftWallPmirror);
%     Pold(Data.RightWallP)=Pold(Data.RightWallPmirror);
        Pressure(Data.TopWallP)=Pold(Data.TopWallPmirror);
    Pressure(Data.BottomWallP)=Pold(Data.BottomWallPmirror);
    Pressure(Data.LeftWallP)=Pold(Data.LeftWallPmirror);
    Pressure(Data.RightWallP)=Pold(Data.RightWallPmirror);

    for i = (1:xSize)
        for j = (1:ySize)
            if IsCenterP(j,i)==true %checks if node is central node
                 Pressure(j,i) = (1-SOR).*Pold(j,i)+SOR.*(BetaSquared*(Pressure(j-1,i)+Pold(j+1,i))+(Pressure(j,i-1)+Pold(j,i+1))-C(j,i)*dx^2)/Divisor;
%                  Pressure(j,i) = (1-SOR).*Pold(j,i)+SOR.*(Pold(j,i+1)+Pressure(j,i-1)+Beta^2.*(Pold(j+1,i)+Pressure(j-1,i))-dx^2.*ConstantMat(j,i))./(2*(1+Beta^2));
           

                
            end
        end
    end
    
    Error2 = norm(Pressure-Pold,'fro'); %Calculate norm 2 error
    if Iterations ==50000
        Stop=1;
    end
    Pold=Pressure;
    Iterations = Iterations+1;
end
Done=1;
end

