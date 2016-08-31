function powerChannel = PowerPerChannel(wlan,numChannels,powMat)
    powMat = db2pow(powMat);
    powerChannel = zeros(size(wlan,2),numChannels);
    for j=1:size(wlan,2)
        for k=1:size(wlan,2)
            for i=1:numChannels
                if j~=k && wlan(k).channel == i %&& wlan(k).transmitting
                    powerChannel(j,i) = powerChannel(j,i) + powMat(j,k);
                end
            end
        end
    end
end