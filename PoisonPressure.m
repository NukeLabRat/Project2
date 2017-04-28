function [ P, Iterations] = PoisonPressure( ConstantMat, IsCenter)
%PoisonPressure Pressure solving function
%   Itteratively solves for the pressure field durring each timestep. Gives
%   back the pressure field in a matrix at locations given in NodeX and
%   NodeY. 


P=zeros(yEnd,xEnd);
Pold=P;
Iterations = 0;
Error2 = 1;
SOR=1; %1.7189 is optimal value.
while Error2>1E-8
    for i = 1:xEnd
        for j = 1:yEnd
            if IsCenter(j,i)==true %checks if node is central node
            P(j,i) = (1-SOR).*Pold(j,i)+SOR.*(Pold(j,i+1)+P(j,i-1)+Beta^2.*(Pold(j+1,i)+P(j-1,i))+dx^2.*ConstantMat)./4;
            
            else %For Boundary Nodes
                if i==1
                    P(j,i-1)=Pold(j,i);
                end
                if i==xEnd
                    Pold(j,i+1)=Pold(j,i);
                end
                if j==1
                    P(j-1,i)=Pold(j,i);
                end
                if j==yEnd
                    Pold(j+1,i)=Pold(j,i);
                end
                
                    P(j,i) = (1-SOR).*Pold(j,i)+SOR.*(Pold(j,i+1)+P(j,i-1)+Beta^2.*(Pold(j+1,i)+P(j-1,i))+dx^2.*ConstantMat)./4;
                
            end
            
        end
    end
    
    Error2 = norm(P-Pold); %Calculate norm 2 error
    Pold=P;
    Iterations = Iterations+1;
end


end

