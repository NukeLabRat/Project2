function [ P ] = PoisonPressure( ConstantMat, NodeX, NodeY, CenterNodes)
%PoisonPressure Pressure solving function
%   Itteratively solves for the pressure field durring each timestep. Gives
%   back the pressure field in a matrix at locations given in NodeX and
%   NodeY. 
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

k = 2;
Error2 = 1;
Error1 = 1;
tic
SOR=1.7189; %1.7189 is optimal value.
while Error2(k-1)>1E-8
    for i = 1:xEnd
        for j = 1:yEnd
 
            if IsCenter(j,i)==true %checks if node is central node
            Psi(j,i,k) = (1-SOR).*Psi(j,i,k-1)+SOR.*(Psi(j,i+1,k-1)+Psi(j,i-1,k)+Beta^2.*(Psi(j+1,i,k-1)+Psi(j-1,i,k)))./(2*(1+Beta^2));
            else
                Psi(j,i,k) = Psi(j,i,k-1); %applies BC condition of Boundary Nodes
            end
            
        end
    end
    Error1(k-1) = norm(Psi(:,:,k)-Psi(:,:,k-1),1); %Calculate norm 1 error
    Error2(k) = norm(Psi(:,:,k)-Psi(:,:,k-1)); %Calculate norm 2 error
    ErrorInf(k-1) = norm(Psi(:,:,k)-Psi(:,:,k-1),Inf); %Calculate norm 3 error
    k = k+1;
end
toc
end

