function T=ThomasMat(Mat,d)
% Thomas solver for tridagonal matrix - expects matrix input
%
% a is subdiagonal vector
% b is the diagonal vector
% c is the superdiagonal vector
% a and c must be the same length as b, just pad them both appropriately.
% The values a(1) and b(N) are never called and thus can be anything;

N=length(d);
b=diag(Mat);
a(1) = 1; %padding vector
a(2:N) = diag(Mat,-1);
c = diag(Mat,1);
c(N) = 1; %padding vector

for Cnodes = 2:N; %Central Nodes
b(Cnodes) = b(Cnodes)-a(Cnodes)./b(Cnodes-1).*c(Cnodes-1);
d(Cnodes) = d(Cnodes)-a(Cnodes)./b(Cnodes-1).*d(Cnodes-1);
end

T(N) = d(N)./b(N);
for Nodes = N-1:-1:1; %Backsubing in reverse
T(Nodes) = (d(Nodes)-c(Nodes).*T(Nodes+1))./b(Nodes);
end
end