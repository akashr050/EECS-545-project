% author: Yutong Wang

function A = grid_graph( n_squared )
% Makes an sqrt(n_squared)-by-sqrt(n_squared) grid matrix
% input: n_squared an square integer
% output: a graph adjacency matrix without ones on the diagonal

n = sqrt(n_squared);
if mod(n,1)~= 0
    error('grid_graph requires a square number of features');
end
A = zeros(n^2, n^2);

toIntRepr = @(t) (t(1)-1)*n+t(2);

for a = 1:n
    for b = 1:n
        s = toIntRepr([a,b]);
        if a > 1
            t = toIntRepr([a-1,b]);
            A(s,t) = 1;
            A(t,s) = 1;
        end
        if a < n
            t = toIntRepr([a+1,b]);
            A(s,t) = 1;
            A(t,s) = 1;
        end
        if b > 1
            t = toIntRepr([a,b-1]);
            A(s,t) = 1;
            A(t,s) = 1;
        end
        if b < n
            t = toIntRepr([a,b+1]);
            A(s,t) = 1;
            A(t,s) = 1;
        end
    end
end
end
