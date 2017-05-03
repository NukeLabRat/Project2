function [ Pressure, Iterations] = PoisonPressureSLOR3(Data)
%PoisonPressure Pressure solving function
%   Itteratively solves for the pressure field durring each timestep. Gives
%   back the pressure field in a matrix at locations given in NodeX and
%   NodeY.


Iterations = 0;
Error2 = 1;
Beta=Data.dx/Data.dy;
BetaSquared=Beta^2;
Pold=Data.P0;
Pressure=Pold;
C=Data.ConstantMat.*Data.dx.^2; %Constant given by ustar and vstar
Divisor=2.*BetaSquared+2;
b = zeros(Data.xSize,1);
PressureMat = Data.P0;
Pold=Data.P0;
SOR=1.1;
while Error2>1E-8
    
    Pressure(Data.TopWallP)=Pold(Data.TopWallPmirror);
    Pressure(Data.BottomWallP)=Pold(Data.BottomWallPmirror);
    Pressure(Data.LeftWallP)=Pold(Data.LeftWallPmirror);
    Pressure(Data.RightWallP)=Pold(Data.RightWallPmirror);
    
    for j = 1:Data.ySize
        for i = 1:Data.xSize
            if Data.IsCenterP(j,i)==true
                PressureMat(i,i) = Divisor;
                PressureMat(i,i+1) = -SOR;
                PressureMat(i,i-1) = -SOR;
                b(i) = Divisor.*(1-SOR).*Pold(j,i)+SOR.*(Beta^2.*(Pold(j+1,i)+Pressure(j-1,i))+C(j,i));
            else
                PressureMat(i,i) = 1;
                b(i) = Pold(j,i);
            end
            
        end
        Pressure(j,:) = ThomasMat(PressureMat,b);
%         Pressure(j,:) = PressureMat\b;         
    end
    
    
   
    Error2 = norm(Pressure(Data.Center)-Pold(Data.Center),'fro'); %Calculate norm 2 error
    if Iterations ==5000
        Stop=1; %Place to put breakpoint when debugging.
    end
    Pold=Pressure;
    Iterations = Iterations+1;
end
% Done=1;
end

