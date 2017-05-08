function [ Pressure, Iterations] = PoisonPressureSLOR4(Data)
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
SOR=Data.SOR;
while Error2>Data.Error
    
%     Pressure(Data.TopWallP)=Pold(Data.TopWallPmirror);
%     Pressure(Data.BottomWallP)=Pold(Data.BottomWallPmirror);
%     Pressure(Data.LeftWallP)=Pold(Data.LeftWallPmirror);
%     Pressure(Data.RightWallP)=Pold(Data.RightWallPmirror);
            Pold(:,1) = Pold(:,2);
        Pold(:,end) = Pold(:,end-1);
        Pold(1,:) = Pold(2,:);
        Pold(end,:) = Pold(end-1,:);
    
    for j = 1:Data.ySize
        for i = 1:Data.xSize
            if Data.IsCenterP(j,i)==true
                PressureMat(i,i) = Divisor;
                PressureMat(i,i+1) = -SOR;
                PressureMat(i,i-1) = -SOR;
                b(i) = Divisor.*(1-SOR).*Pold(j,i)+SOR.*(Beta^2.*(Pold(j+1,i)+Pressure(j-1,i))-C(j,i));
            else
                PressureMat(i,i) = 1;
                b(i) = Pold(j,i);
            end
            
        end
        Pressure(j,:) = ThomasMat(PressureMat,b);
%         Pressure(j,:) = PressureMat\b;         
    end
    
    
   
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
% Done=1;
end

