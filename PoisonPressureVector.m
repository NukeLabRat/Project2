function [ Pressure, Iterations] = PoisonPressureVector(Data)
%PoisonPressure Pressure solving function
%   Itteratively solves for the pressure field durring each timestep. Gives
%   back the pressure field in a matrix at locations given in NodeX and
%   NodeY.


Iterations = 0;
Error2 = 1;
SOR=1; %1.7189 is optimal value.
Beta=Data.dx/Data.dy;
% [ySize, xSize] = size(IsCenterP);
% ConstantMat(isnan(ConstantMat))=0;
Pold=Data.P0;
Pressure=Pold;
while Error2>1E-8
    
    Pressure(Data.BottomWallP)=Pold(Data.BottomWallPmirror);
    Pressure(Data.LeftWallP)=Pold(Data.LeftWallPmirror);
    Pressure(Data.RightWallP)=Pold(Data.RightWallPmirror);
    Pressure(Data.TopWallP)=Pold(Data.TopWallPmirror);
    
    Pressure(Data.Center) = (1-SOR).*Pold(Data.Center)+SOR.*(Pold(Data.Iplus)+Pold(Data.Iminus)+Beta.^2.*(Pold(Data.Jplus)+Pold(Data.Jminus))-Data.dx.^2.*Data.ConstantMat(Data.Center))./(2.*(1+Beta.^2));
    
    
    Error2 = norm(Pressure-Pold); %Calculate norm 2 error
    if Iterations ==50000
        Stop=1; %Place to put breakpoint when debugging.
    end
    Pold=Pressure;
    Iterations = Iterations+1;
end
% Done=1;
end

