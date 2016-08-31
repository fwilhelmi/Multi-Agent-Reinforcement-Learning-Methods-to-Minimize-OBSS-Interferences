function selected_action = ChooseAction(Qval,curr_state,actions_ch,actions_cca,actions_tpc)
%choose_action: returns the best possible action given the current state
%   Inputs:
%       - Qval :
%       - curr_state :
%       - actions_ch :
%       - actions_cca :
%       - actions_tpc :
%       - e : Epsilon-greedy approach for exploration
%   Returns "selected_action", which contains the chosen channel, CCA and TPC
    indexes=[];

    [val,~] = max(Qval);
    % Check if there is more than one occurrence in order to select a value randomly
    if sum(Qval(:)==val)>1
        for i=1:size(Qval,2)
            if Qval(i) == val, indexes = [indexes i]; end
        end
        index = randsample(indexes,1);
    else
        [~,index] = max(Qval);
    end

    [a,b,c] = val2indexes(index,size(actions_ch,2),size(actions_cca,2),size(actions_tpc,2));
    selected_action = [a b c];
end