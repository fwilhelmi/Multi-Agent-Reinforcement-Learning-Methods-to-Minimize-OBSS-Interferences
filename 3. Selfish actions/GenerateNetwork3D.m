function wlan = GenerateNetwork3D(N_WLANs, NumChannels, B)
    % DrawNetwork3D  Calculate interferences on WLANs.
    %   Inputs:
    %       * N_WLANs: number of WLANs on the studied environment
    %       * NumChannels: number of available channels
    %       * B: bandwidth available per WLAN (Hz)
    %   Output:
    %       * wlan: object containing the information of each wlan drawn
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OTHER CONFIGURABLE PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    min_STA = 1; %minimum number of associated STAs
    max_STA = 4; %maximum number of associated STAs
    % Number of associated STAs on each WLAN
    N_STAs = round((max_STA - min_STA).*rand(1,N_WLANs) + min_STA);
    %disp('STAs per AP');
    %disp(N_STAs)
    printSTAs = 0;
    %colors = {[0 0 1],[0 0 .5],[0 0.5 0],[0 1 0],[.5 0 0],[1 0 0],[0.2 0.2 1],[1 .2 .2],[.2 1 .2],[.3 .8 .2]};

%     % Change the seed for the random generators 'rand' and 'randn'
%     seed = 4;
%     rand('seed', seed); 
%     randn('seed', seed); 
    % - rand: random elements uniformly distributed on the interval (0, 1)
    % - randn: normally distributed random elements
    % Dimensions of the 3D map
    MaxX=10;
    MaxY=5; 
    MaxZ=10;
    % Maximum range for a STA
    MaxRangeX = 5;
    MaxRangeY = 5;
    MaxRangeZ = 1;
    % AP density
    disp('Density of APs');
    disp(N_WLANs/(MaxX*MaxY*MaxZ));
 
    %% Locate elements on the map randomly
    for j=1:N_WLANs    
        % Assign Tx Power and CCA on the WLAN
        wlan(j).PTdBm = 20;
        wlan(j).CCA = -82;
        % Assign channel to the AP randomly
        wlan(j).channel = ceil(NumChannels*rand());
        % Assign location to the AP on the 3D map
        wlan(j).x = MaxX*rand();
        wlan(j).y = MaxY*rand();
        wlan(j).z = MaxZ*rand();  
        % Build arrays of locations for each AP
        x(j)=wlan(j).x;
        y(j)=wlan(j).y;
        z(j)=wlan(j).z;
        % Assign a STA to each AP for throughput analysis
        if(rand() < 0.5), xc = MaxRangeX.*rand();   %dnode*rand(); %what xc represents ? B: Is just an auxiliary variable to fix the position of the node around the AP, see below
        else xc = -MaxRangeX.*rand();
        end
        if(rand() < 0.5), yc = MaxRangeY.*rand();
        else yc = -MaxRangeY.*rand();
        end
        if(rand() < 0.5), zc = MaxRangeZ.*rand();
        else zc = -MaxRangeZ.*rand();
        end
        wlan(j).xn = min(abs(wlan(j).x+xc), MaxX);  
        wlan(j).yn = min(abs(wlan(j).y+yc), MaxY);
        wlan(j).zn = min(abs(wlan(j).z+zc), MaxZ);
        xn(j)=wlan(j).xn; %what is xn(j) B: the "x" position of node j
        yn(j)=wlan(j).yn;
        zn(j)=wlan(j).zn;        
        % Assign location of STAs associated to each AP
        %  - Example (access to coordinates of each STA):
        %       wlan(1).stas{1}(1) - access to "x" of STA 1 at WLAN 1
        %       wlan(2).stas{4}(3) - access to "z" of STA 4 at WLAN 2
        wlan(j).stas = {};
        maxD = 0;
        maxD_index = 0;
        for k=1:N_STAs(j)
            d = [MaxRangeX*rand() MaxRangeY*rand() MaxRangeZ*rand()];
            if(rand() < 0.5), xs = d(1);
            else xs = -d(1);
            end    
            if(rand() < 0.5), ys = d(2);
            else ys = -d(2);
            end     
            if(rand() < 0.5), zs = d(3);
            else zs = -d(3);
            end
            xs = min(abs(wlan(j).x+xs), MaxX);  
            ys = min(abs(wlan(j).y+ys), MaxY);
            zs = min(abs(wlan(j).z+zs), MaxZ);
            wlan(j).stas = [wlan(j).stas; [xs ys zs]]; 
            %sta(j,k)
            d_AP_STA = sqrt(([x(j) y(j) z(j)] - [xs ys zs]).^2);
            if d_AP_STA > maxD
                maxD = d_AP_STA;
                maxD_index = k;
            end
        end
        % Maximum distance from an AP to its associated STAs
        wlan(j).maxDistance = maxD;
        wlan(j).maxDistanceSta = maxD_index;
        wlan(j).BW = B; 
    end
    disp('Channels selected per WLAN')
    for i=1:N_WLANs
        channels(i) = wlan(i).channel;       
    end
    for i=1:NumChannels
        sumChannels(i) = sum(channels==i);
    end
    disp(channels)
    disp('Times a channel is occupied')
    disp(sumChannels)
    %% Plot map of APs and STAs
   if printSTAs == 1
        for j=1:size(wlan,2)
            x(j)=wlan(j).x;
            y(j)=wlan(j).y;
            z(j)=wlan(j).z;
        end
        figure
        axes;
        set(gca,'fontsize',12);
        labels = num2str((1:size(y' ))','%d');  
        for i=1:N_WLANs
            color = [rand() rand() rand()];
            scatter3(wlan(i).x,wlan(i).y,wlan(i).z,50,color,'filled');
            hold on;
            scatter3(wlan(i).xn,wlan(i).yn,wlan(i).zn,30,color,'filled');
            line([wlan(i).x,wlan(i).xn],[wlan(i).y,wlan(i).yn],[wlan(i).z,wlan(i).zn],'Color',color);
%             for j=1:N_STAs(i)
%                 scatter3(wlan(i).stas{j}(1),wlan(i).stas{j}(2),wlan(i).stas{j}(3),5,[0 0 0],'filled');
%                 line([wlan(i).x,wlan(i).stas{j}(1)],[wlan(i).y,wlan(i).stas{j}(2)],[wlan(i).z,wlan(i).stas{j}(3)],'Color',color);
%                 hold on;
%             end
        end
        axis([0 MaxX 0 MaxY 0 MaxZ])
        text(x,y,z,labels,'horizontal','left','vertical','bottom') 
        xlabel('x [meters]','fontsize',12);
        ylabel('y [meters]','fontsize',12);
        zlabel('z [meters]','fontsize',12);
    end