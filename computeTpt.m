function tpt = computeTpt(wlan,powMat,time_slots,noise)
    N_WLANs = size(wlan,2);
    % COMPUTE SNR AND TPT FOR EACH WLAN ACCORDING TO BACKOFF SIMULATION IN ORDER TO PROVIDE A JOINT REWARD
    avg_transmission_rate = zeros(N_WLANs,time_slots);
    WLAN_is_transmitting = zeros(1,N_WLANs); % We store here the number of times an AP transmits
    t = 1;
    while t < time_slots            
        for i=1:N_WLANs, wlan(i).transmitting = 0; end
        BO_values = rand(1,N_WLANs);
        for i=1:N_WLANs
            [a,b]=min(BO_values); % a the value, b its position in the vector
            BO_values(b)=Inf; % to remove it
            wlan(b).transmitting = 1;
            interferences = Interferences(wlan, powMat); % mW
            %BO_values_stas = rand(1,wlan(b).STAs);
            if(pow2db((interferences(b)+db2pow(noise))) < wlan(b).CCA), WLAN_is_transmitting(b)=WLAN_is_transmitting(b)+1;
            else wlan(b).transmitting = 0; end
        end
        interferences = Interferences(wlan, powMat); %dBm
        snr_aux = zeros(1,N_WLANs);
        for i=1:N_WLANs
            snr_aux(i) = powMat(i,i) - pow2db((interferences(i)+db2pow(noise))); % dBm
            avg_transmission_rate(i,t) = Capacity(wlan(i).BW, db2pow(snr_aux(i))); % bps
        end
        t = t+1;
    end
    tpt = zeros(1,N_WLANs);
    for i=1:N_WLANs
        tpt(i) = mean(avg_transmission_rate(i,:))*(WLAN_is_transmitting(i)/time_slots); % bps
    end
end