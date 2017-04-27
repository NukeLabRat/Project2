function [ Edges ] = MatEdges( A )
%MatEdges Returns vector of of "edges" of matrix
%   Given a matrix, MatEdges returns a vector of all elements that are on a
%   "boundary". Used for applying boundary conditions;

[y, x] = size(A);
Edge1 = sub2ind([y x],1:y,ones(1,y));
Edge2 = sub2ind([y x],1:y,x*ones(1,y));
Edge3 = sub2ind([y x],ones(1,x),1:x);
Edge4 = sub2ind([y x],y*ones(1,x),1:x);
EdgeIndices = [Edge1 Edge2 Edge3 Edge4];
Edges=unique(EdgeIndices);
end

