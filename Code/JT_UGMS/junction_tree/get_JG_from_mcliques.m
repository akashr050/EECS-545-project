% Author: Yutong Wang

function JG = get_JG_from_mcliques( mcliques )
%get_jg_from_mcliques Generates the junction graph
%  input - mcliques, n x p the rows are 0-1 matrices, n = number of
%  cliques, p = number of features
%  output - JG, junction graph
    JG = mcliques*(mcliques');
    JG = triu(JG,1) + tril(JG,-1);
end

