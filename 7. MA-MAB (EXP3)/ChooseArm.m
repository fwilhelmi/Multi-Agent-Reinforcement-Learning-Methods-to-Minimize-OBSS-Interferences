function arm = ChooseArm(p_per_arm)
%ChooseArm: returns the best possible arm given the current distribution
%   Inputs:
%       - rewards_per_configuration: rewards noticed at each configuration
%       - e: epsilon
%   Returns "arm", which represents the configuration composed by [channel,CCA,TPC]
    indexes=[];
    
    randsample(1:24, 1, true, p_per_arm)
    
    % Exploration approach
    if rand()>e 
        [val,~] = max(rewards_per_configuration);
        % Break ties randomly
        if sum(rewards_per_configuration==val)>1
            if val ~= Inf
                indexes = find(rewards_per_configuration==val);
                arm = randsample(indexes,1);
            else
                arm = randsample(1:size(rewards_per_configuration,2),1);
            end
        % Select arm with maximum reward
        else
            [~,arm] = max(rewards_per_configuration);
        end
    else
        arm = randi([1 size(rewards_per_configuration,2)], 1, 1);
    end
end