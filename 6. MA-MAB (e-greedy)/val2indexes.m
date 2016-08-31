function [i,j,k] = val2indexes(x,a,b,c)
% We can know i,j,k of each states with this (e.g. state x)
%   k = ceil(x/(size(actions_TxPower,2)*size(actions_CCA,2)); 
%   j = ceil(x/(size(actions_channel,2)); 
%   i = mod(x,size(actions_channel,2)+1); -> obtaining 0 means max(actions_channel)
    i = mod(x,a); 
    if i == 0, i = a; end    
    y = mod(x,(a*b));
    j = ceil(y/a);
    if j == 0, j = b; end 
    k = ceil(x/(a*b)); 
    if k > c, k = c; end
end
