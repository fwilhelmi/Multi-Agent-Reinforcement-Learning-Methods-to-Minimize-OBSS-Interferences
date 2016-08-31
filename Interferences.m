function interferences = Interferences(wlan, powMat)
%%Interferences - Returns the interferences power received at each WLAN
%   Input:
%       - wlan: contains information of each WLAN in the map. For instance,
%       wlan(1) corresponds to the first one, so that it has unique
%       parameters (x,y,z,BW,CCA,etc.).
%       - powMat: matrix NxN (N is the number of WLANs) with the power
%       received at each AP in dBm.
%   Output:
%       - intMat: 1xN array (N is the number of WLANs) with the
%       interferences noticed on each AP in mW
% We assume that overlapping channels also create an interference with lower level (20dB/d) 
% 20 dB == 50 dBm
    interferences = zeros(1,size(wlan,2));  
    for i=1:size(wlan,2)
        for j=1:size(wlan,2)
            if i~=j && wlan(j).transmitting == 1
                interferences(i) = interferences(i) + db2pow(powMat(i,j) - db2pow(50*(abs(wlan(i).channel-wlan(j).channel))));
            end
        end
    end
%     intMat = pow2db(intMat);
end