function [ P, Iterations] = PoisonPressure( ConstantMat, CenterNodes, yEnd, xEnd, UstarX, UstarY)
%PoisonPressure Pressure solving function
%   Itteratively solves for the pressure field durring each timestep. Gives
%   back the pressure field in a matrix at locations given in NodeX and
%   NodeY. 

Index=zeros(yEnd-1,xEnd-1);
IsCenter=zeros(yEnd-1,xEnd-1);
    for i = 1:xEnd-1 %This loop determines whether each node is central node or boundary node.
        for j = 1:yEnd-1
            Index(j,i) = sub2ind([yEnd-1 xEnd-1],j,i);
            if sum(ismember(CenterNodes,Index(j,i)))==1
                IsCenter(j,i) = true;
            else
                IsCenter(j,i) = false;
            end
        end
    end
    rng('default')
P=rand(yEnd+1,xEnd+1);
ConstantMat=rand(size(P));
Pold=P;
Iterations = 0;
Error2 = 1;
SOR=1; %1.7189 is optimal value.
Beta=1;
dx=.25;
while Error2>1E-8
    for i = (2:xEnd)
        for j = (2:yEnd)
            if IsCenter(j-1,i-1)==true %checks if node is central node
            P(j,i) = (1-SOR).*Pold(j,i)+SOR.*(Pold(j,i+1)+P(j,i-1)+Beta^2.*(Pold(j+1,i)+P(j-1,i))+dx^2.*ConstantMat(j,i))./(2*(1+Beta^2));
            
            else %For Boundary Nodes
                if i==2
                    P(j,i-1)=Pold(j,i);
                end
                if i==xEnd
                    Pold(j,i+1)=Pold(j,i);
                end
                if j==2
                    P(j-1,i)=Pold(j,i);
                end
                if j==yEnd
                    Pold(j+1,i)=Pold(j,i);
                end
                    P(j,i) = (1-SOR).*Pold(j,i)+SOR.*(Pold(j,i+1)+P(j,i-1)+Beta^2.*(Pold(j+1,i)+P(j-1,i))+dx^2.*ConstantMat(j,i))./(2*(1+Beta^2));
                
            end
            
        end
    end
    
    Error2 = norm(P(2:17,2:25)-Pold(2:17,2:25)); %Calculate norm 2 error
    if Iterations ==5000
        h=1;
    end
    Pold=P;
    Iterations = Iterations+1;
end


end

