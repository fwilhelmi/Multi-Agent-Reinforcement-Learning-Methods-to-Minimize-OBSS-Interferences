function [action, Tpt] = bestAction(wlan, index, parameter, feasible_actions)
    noise = -100;
    rewards = zeros(1,length(feasible_actions));
    curr_state = [wlan(index).channel, wlan(index).CCA, wlan(index).PTdBm];
    curr_tpt = wlan(index).tpt;
    wlan_aux = wlan;
    for i=1:length(feasible_actions)
        if parameter == 1 %Choose better channel
            wlan_aux(index).channel = feasible_actions(i);    
        elseif parameter == 2 %Choose better CCA level
            wlan_aux(index).CCA = feasible_actions(i);    
        elseif parameter == 3 %Choose better TPC level
            wlan_aux(index).PTdBm = feasible_actions(i);  
        end
        powMat = PowerMatrix(wlan_aux); %                    
        interferences = Interferences(wlan_aux, powMat); %dBm
        % Calculate SNR and Tpt experienced by the WLAN
        wlan_aux(index).snr = powMat(index,index) - pow2db((interferences(index)+db2pow(noise))); % dBm
        wlan_aux(index).tpt = Capacity(wlan_aux(index).BW, db2pow(wlan_aux(index).snr)); % bps
        % Save obtained tpt on rewards
        rewards(i) = wlan_aux(index).tpt;
    end
    [bestR indexR] = max(rewards);
%     action = feasible_actions(indexR);
%     Tpt = bestR;
    if bestR < curr_tpt
        action = curr_state(parameter);
        Tpt = curr_tpt;
    else
        action = feasible_actions(indexR);
        Tpt = bestR;
    end   
end

