function DrawNetwork3D(wlan)
%%DrawNetwork 3D - Plots a 3D of the network and highlights the
%%interfering WLANs
%   Input:
%       - wlan: contains information of each WLAN in the map. For instance,
%       wlan(1) corresponds to the first one, so that it has unique
%       parameters (x,y,z,BW,CCA,etc.).
%       - powMat: matrix NxN (N is the number of WLANs) with the power
%       received at each AP in dBm.
    % Size of the map
    MaxX=8;
    MaxY=6; 
    MaxZ=15;
    for j=1:size(wlan,2)
        x(j)=wlan(j).x;
        y(j)=wlan(j).y;
        z(j)=wlan(j).z;
    end
    figure
    axes;
    set(gca,'fontsize',12);
    labels = num2str((1:size(y' ))','%d');    
    for i=1:size(wlan,2)
        scatter3(wlan(i).x,wlan(i).y,wlan(i).z,50,[0 0 0],'filled');
        hold on;   
        %axis([-10 MaxX -10 MaxY -10 MaxZ]);
        %labels_sta = {};
        for j=1:size(wlan(i).stas,1)
            xn=wlan(i).stas{j}(1);
            yn=wlan(i).stas{j}(2);
            zn=wlan(i).stas{j}(3);
            scatter3(xn,yn,zn,20,[0 0 1],'filled');
            hold on;
            line([wlan(i).x,wlan(i).stas{j}(1)],[wlan(i).y,wlan(i).stas{j}(2)],[wlan(i).z,wlan(i).stas{j}(3)],'Color',[0.4,0.4,1.0],'LineStyle',':');
%             labels_sta = [labels_sta;['STA ' num2str(i) '.' num2str(j)]];        
%             text(xn,yn,zn,labels_sta,'horizontal','left','vertical','bottom')
        end
    end
    text(x,y,z,labels,'horizontal','left','vertical','bottom') 
    xlabel('x [meters]','fontsize',12);
    ylabel('y [meters]','fontsize',12);
    zlabel('z [meters]','fontsize',12);
end