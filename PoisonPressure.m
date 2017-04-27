function [ P, Iterations] = PoisonPressure( ConstantMat, NodeX, NodeY, CenterNodes)
%PoisonPressure Pressure solving function
%   Itteratively solves for the pressure field durring each timestep. Gives
%   back the pressure field in a matrix at locations given in NodeX and
%   NodeY. 


P=zeros(yEnd,xEnd);
Index=zeros(yEnd,xEnd);
IsCenter=zeros(yEnd,xEnd);
    for i = 1:xEnd %This loop determines whether each node is central node or boundary node.
        for j = 1:yEnd
            Index(j,i) = sub2ind([yEnd xEnd],j,i);
            if sum(ismember(CenterNodes,Index(j,i)))==1
            IsCenter(j,i) = true;
            else
                IsCenter(j,i) = false;
            end
        end
    end

Iterations = 0;
Error2 = 1;
SOR=1; %1.7189 is optimal value.
while Error2>1E-8
    for i = 1:xEnd
        for j = 1:yEnd
            if IsCenter(j,i)==true %checks if node is central node
            P(j,i) = (1-SOR).*Pold(j,i)+SOR.*(Pold(j,i+1)+P(j,i-1)+Beta^2.*(Pold(j+1,i)+P(j-1,i))+dx^2.*ConstantMat)./4;
            else
                P(j,i) = Pold(j,i); %applies BC condition of Boundary Nodes
            end
            
        end
    end
    
    Error2 = norm((P-Pold); %Calculate norm 2 error
    Pold=P;
    Iterations = Iterations+1;
end

P=P
end

