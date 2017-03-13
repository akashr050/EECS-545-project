function Sij = SeparatorSearchSpaceNew(nei,nej,i,j,H,pcType,eta)


if (length(nei) < length(nej))
    ne = [i j mysetdiff(nei,j)];
else
    ne = [i j mysetdiff(nej,i)];
end

if (length(ne)-1 < eta)
    Sij = [];
    Sij2 = [];
    disp('sadfdsfdsfsfds')
    return;
end

Sij = setdiff(ne,[i j]);

end