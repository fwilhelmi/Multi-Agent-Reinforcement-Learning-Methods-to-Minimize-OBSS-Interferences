function x = indexes2val(i,j,k,a,b)
% We can know the state x from i,j,k (value for each variable)
    x = i + (j-1)*a + (k-1)*a*b;
end