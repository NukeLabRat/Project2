function [ Pressure, Iterations] = PoisonPressure( ConstantMat, IsCenterP, P0)
%PoisonPressure Pressure solving function
%   Itteratively solves for the pressure field durring each timestep. Gives
%   back the pressure field in a matrix at locations given in NodeX and
%   NodeY. 




Iterations = 0;
Error2 = 1;
SOR=1; %1.7189 is optimal value.
Beta=1;
dx=.25;
[ySize, xSize] = size(IsCenterP);
Pzero=zeros(ySize+2,xSize+2);
Pzero(2:ySize+1,2:xSize+1)=P0;
Pold=Pzero;
P=Pold;
while Error2>1E-8
    for i = (2:xSize+1)
        for j = (2:ySize+1)
            if IsCenterP(j-1,i-1)==true %checks if node is central node
            P(j,i) = (1-SOR).*Pold(j,i)+SOR.*(Pold(j,i+1)+P(j,i-1)+Beta^2.*(Pold(j+1,i)+P(j-1,i))+dx^2.*ConstantMat(j,i))./(2*(1+Beta^2));
            
            else %For Boundary Nodes
                if i==2
                    P(j,i-1)=Pold(j,i);
                end
                if i==xSize+1
                    Pold(j,i+1)=Pold(j,i);
                end
                if j==2
                    P(j-1,i)=Pold(j,i);
                end
                if j==ySize+1
                    Pold(j+1,i)=Pold(j,i);
                end
                    P(j,i) = (1-SOR).*Pold(j,i)+SOR.*(Pold(j,i+1)+P(j,i-1)+Beta^2.*(Pold(j+1,i)+P(j-1,i))+dx^2.*ConstantMat(j,i))./(2*(1+Beta^2));
                
            end
            
        end
    end
    
    Error2 = norm(P(2:ySize+1,2:xSize+1)-Pold(2:ySize+1,2:xSize+1)); %Calculate norm 2 error
    if Iterations ==5000
        h=1;
    end
    Pold=P;
    Iterations = Iterations+1;
end

Pressure=P(2:ySize+1,2:xSize+1);
end

