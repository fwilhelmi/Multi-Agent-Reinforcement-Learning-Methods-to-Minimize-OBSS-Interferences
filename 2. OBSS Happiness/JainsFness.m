function [ fness ] = JainsFness(C)
%jains_fairness calculates the Jain's fairness measure given the input:
%   input (C): array of capacities that each WLAN experiences
%   output (fness): Jain's fairness measure [0,1]
    fness = sum(C)^2 ./ (size(C,2)*sum(C.^2));
%     disp('Jains fairness coefficient:');
%     disp(fness);
end