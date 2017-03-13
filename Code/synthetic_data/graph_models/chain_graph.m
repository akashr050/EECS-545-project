% Author: Yutong Wang

function A = chain_graph( n )
%AR(1) model from http://pages.stat.wisc.edu/~myuan/papers/graph.final.pdf

A = eye(n);
for i =1:n-1
    A(i+1,i) = 0.5;
    A(i,i+1)=0.5;
end
end

