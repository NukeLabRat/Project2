function [ Pressure, Iterations] = TestPoisonPressure( ConstantMat, IsCenterP, P0)
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
ConstantMat=zeros(ySize,xSize);
P0=ConstantMat;
Pold=rand(ySize,xSize);
Pold(MatEdges(Pold))=1;
P=Pold;
while Error2>1E-8
    for i = (1:xSize)
        for j = (1:ySize)
            if IsCenterP(j,i)==true %checks if node is central node
                if i==1
                    P(j,i-1)=Pold(j,i);
                end
                if i==xSize
                    Pold(j,i+1)=Pold(j,i);
                end
                if j==1
                    P(j-1,i)=Pold(j,i);
                end
                if j==ySize
                    Pold(j+1,i)=Pold(j,i);
                end
            P(j,i) = (1-SOR).*Pold(j,i)+SOR.*(Pold(j,i+1)+P(j,i-1)+Beta^2.*(Pold(j+1,i)+P(j-1,i))+dx^2.*ConstantMat(j,i))./(2*(1+Beta^2));
            
            else %For Boundary Nodes

            P(j,i) = Pold(j,i);
                
            end
            
        end
    end
    
    Error2 = norm(P-Pold); %Calculate norm 2 error
    if Iterations ==5000
        h=1;
    end
    Pold=P;
    Iterations = Iterations+1;
end

Pressure=P;
end

