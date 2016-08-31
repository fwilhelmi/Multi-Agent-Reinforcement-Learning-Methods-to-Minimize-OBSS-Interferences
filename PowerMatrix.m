function powMat = PowerMatrix(wlan)
%%PowerMatrix - Returns the power received by each AP from all the others
%   Input:
%       - wlan: contains information of each WLAN in the map. For instance,
%       wlan(1) corresponds to the first one, so that it has unique
%       parameters (x,y,z,BW,CCA,etc.).
%   Output:
%       - powMat: matrix NxN (N is the number of WLANs) with the power
%       received at each AP in dBm.
    N_WLANs = size(wlan,2); %Number of WLANs (obtained from the input)
    PLd1=5; %Path-loss factor
    shadowing = 9.5; %Shadowing factor
    obstacles = 30; %Obstacles factor
    shadowingmatrix = shadowing*randn(N_WLANs); %Shadowing affecting each WLAN
    obstaclesmatrix = obstacles*rand(N_WLANs); %Obstacles affecting each WLAN
    % Compute the received power on all the APs from all the others
    for i=1:N_WLANs
        for j=1:N_WLANs
            if(i~=j)
                %distance between APs of interest
                d_AP_AP = sqrt((wlan(i).x-wlan(j).x)^2 + (wlan(i).y-wlan(j).y)^2 + (wlan(i).z-wlan(j).z)^2); 
                % Propagation model
                alfa = 4.4;
                %PL = PLd1 + 10*alfa*log10(d) + shawdowing*randn(1) + (d/10).*obstacles.*rand(1);
                PL_AP = PLd1 + 10*alfa*log10(d_AP_AP) + shadowingmatrix(i,j) + (d_AP_AP/10).*obstaclesmatrix(i,j);
                powMat(i,j) = wlan(j).PTdBm - PL_AP;        
            else
                %Calculate Power received at the STA associated to the AP
                d_AP_STA = sqrt((wlan(i).x-wlan(j).xn)^2 + (wlan(i).y-wlan(j).yn)^2 + (wlan(i).z-wlan(j).zn)^2); 
                % Propagation model
                alfa = 4.4;
                %PL = PLd1 + 10*alfa*log10(d) + shawdowing*randn(1) + (d/10).*obstacles.*rand(1);
                PL_AP = PLd1 + 10*alfa*log10(d_AP_STA) + shadowingmatrix(i,j) + (d_AP_STA/10).*obstaclesmatrix(i,j);
                powMat(i,j) = wlan(i).PTdBm - PL_AP;
            end
        end
    end
%     disp('Received Power at each TX ')
%     disp(powMat)   
end