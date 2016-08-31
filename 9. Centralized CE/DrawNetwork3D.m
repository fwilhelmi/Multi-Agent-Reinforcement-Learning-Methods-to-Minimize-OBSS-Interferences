function [] = DrawNetwork3D(wlan,powMat)
    MaxX=10;
    MaxY=6; 
    MaxZ=25;
    for j=1:size(wlan,2)
        x(j)=wlan(j).x;
        y(j)=wlan(j).y;
        z(j)=wlan(j).z;
        xn(j)=wlan(j).xn; %what is xn(j) B: the "x" position of node j
        yn(j)=wlan(j).yn;
        zn(j)=wlan(j).zn;
    end
    figure
    axes;
    set(gca,'fontsize',12);
    scatter3(x,y,z,50,[0 0 0],'filled');
    hold on;        
    labels = num2str((1:size(y' ))','%d');    
    text(x,y,z,labels,'horizontal','left','vertical','bottom') 
    xlabel('x [meters]','fontsize',12);
    ylabel('y [meters]','fontsize',12);
    zlabel('z [meters]','fontsize',12);
    %axis([-10 MaxX -10 MaxY -10 MaxZ]);
    hold on; 
    scatter3(xn,yn,zn,20,[0 0 1],'filled');
    hold on;
    labels_sta = {};
    for i=1:size(wlan,2)
        labels_sta = [labels_sta;['STA ' num2str(i) '.' num2str(1)]];
%         for j=1:size(wlan.STAs,2)
%             labels_sta(i,j) = ['STA ' num2str(i.j)];
%         end
    end
    text(xn,yn,zn,labels_sta,'horizontal','left','vertical','bottom')
    for j=1:size(wlan,2)
        line([wlan(j).x,wlan(j).xn],[wlan(j).y,wlan(j).yn],[wlan(j).z,wlan(j).zn],'Color',[.6 .6 .9]);
        hold on;
        for k=1:size(wlan,2)
            if(powMat(j,k) >= wlan(j).CCA) && j~=k
                %line([wlan(k).x,wlan(j).x],[wlan(k).y,wlan(j).y],[wlan(k).z,wlan(j).z],'Color',[1 0 0]);
            end
        end
    end    
end