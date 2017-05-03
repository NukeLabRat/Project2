function [ Pressure, Iterations] = PoisonPressureVector(Data)
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
C=Data.ConstantMat(Data.Center).*Data.dx.^2; %Constant given by ustar and vstar
Divisor=2.*BetaSquared+2;

while Error2>1E-2
    
    Pressure(Data.TopWallP)=Pold(Data.TopWallPmirror);
    Pressure(Data.BottomWallP)=Pold(Data.BottomWallPmirror);
    Pressure(Data.LeftWallP)=Pold(Data.LeftWallPmirror);
    Pressure(Data.RightWallP)=Pold(Data.RightWallPmirror);
    for i=1:length(Data.Center)
    Pressure(Data.Center(i))= (BetaSquared.*(Pold(Data.Jplus(i))+Pold(Data.Jminus(i)))+Pold(Data.Iplus(i))+Pold(Data.Iminus(i))-C(i))./(Divisor);
    end
    Error2 = norm(Pressure(Data.Center)-Pold(Data.Center),'fro'); %Calculate norm 2 error
    if Iterations ==50000
        Stop=1; %Place to put breakpoint when debugging.
    end
    Pold=Pressure;
    Iterations = Iterations+1;
end
% Done=1;
end

