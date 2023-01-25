%% Best fit method -
% Series of nested for-loops to determine minimum difference orientation
% with coarse pass

function [phi1,theta,phi2]=bestfit1(centcoord)
%Called by CrysChemMain
%For a set of topography coordinates (centcoord) written by topo.m, 
%returns the Euler angle set that minimises angular difference (calls 
%function AngularDiff.m) between orthogonal 5-vector set and centcoord.

if min(size(centcoord))<2 % Topography not indexable when only 1 centcoord 
    % peak is returned by peakfinder function
    [phi1]=NaN;
    [theta]=NaN;
    [phi2]=NaN;
else
    xyz1 = [1 0 0; 0 1 0; 0 0 1;-1 0 0; 0 -1 0]; % Points to rotate
    diffArray = zeros(1,5);
    minDiff = 360;
    minOrient = zeros(1,3);
    centcoord1=centcoord;
    
    centcoord1(:,2) = 90-centcoord(:,2);
    centcoord1 = deg2rad(centcoord1);
    len=length(centcoord);
    centcoord3=ones(len,1);
    centcoord1 = [centcoord1,centcoord3];
    
    % Limits Phi1 1 - 360, theta 1 - 90, Phi2 1 - 360
    % Coarse pass 6 degree increments
    for indexB1=1:60 %Phi 1
        for indexB2=1:15 %Theta
            for indexB3=1:60 %Phi2
                    RBz1=rotz(indexB1*6); % Scale factor to acheive limits 360 degrees
                    RBy=roty(indexB2*6); % " 90 degrees
                    RBz2=rotz(indexB3*6);
                    xyz2 = xyz1*RBz1;
                    xyz3 = xyz2*RBy;
                    xyz4 = xyz3*RBz2;
                    [az,el,r]=cart2sph(xyz4(:,1),xyz4(:,2),xyz4(:,3));
                    xyz4 = horzcat(az,el,r);
                    xyz4(:,1) = rad2deg(xyz4(:,1));
                    centcoord4 = centcoord1;
                    len=size(centcoord4);
    
                    for indexB4=1:length(xyz4)
                        if xyz4(indexB4,1)<0 
                           xyz4(indexB4,1) = xyz4(indexB4,1)+360; %Correcting for circular coordinate convention
                        end
                    end
                    xyz4(:,1) = deg2rad(xyz4(:,1));
                    diff1 = zeros(1,length(centcoord4));
                    for indexB5=1:len(1,1)
                        for indexB6=1:length(xyz4)
                            diffArray(1,indexB6) = AngularDiff(centcoord4(indexB5,:),xyz4(indexB6,:));
                        end
                        diff1(1,indexB5) = min(diffArray);
                    end
                    diffsum = sum(diff1);
                    if diffsum<minDiff
                        minDiff = diffsum;
                        minOrient(1,1) = indexB1*6;
                        minOrient(1,2) = indexB2*6;
                        minOrient(1,3) = indexB3*6;
                    end
            end
        end
    end
    minOrient1=minOrient;
    
    % Fine pass +\- 6 degress from min-orient
    for indexB7=minOrient1(1,1)-6:minOrient1(1,1)+6 %Phi 1
        for indexB8=minOrient1(1,2)-6:minOrient1(1,2)+6 %Theta
            for indexB9=minOrient1(1,3)-6:minOrient1(1,3)+6 %Phi2
                    RBz1=rotz(indexB7);
                    RBy=roty(indexB8);
                    RBz2=rotz(indexB9);
                    xyz2 = xyz1*RBz1;
                    xyz3 = xyz2*RBy;
                    xyz4 = xyz3*RBz2;
                    [az,el,r]=cart2sph(xyz4(:,1),xyz4(:,2),xyz4(:,3));
                    xyz4 = horzcat(az,el,r);
                    xyz4(:,1) = rad2deg(xyz4(:,1));
                    centcoord4 = centcoord1;
                    len=size(centcoord4);
    
                    for indexB4=1:length(xyz4)
                        if xyz4(indexB4,1)<0 
                           xyz4(indexB4,1) = xyz4(indexB4,1)+360; %Correcting for circular coordinate convention
                        end
                    end
                    xyz4(:,1) = deg2rad(xyz4(:,1));
                    diff1 = zeros(1,length(centcoord4));
                    for indexB5=1:len(1,1)
                        for indexB6=1:length(xyz4)
                            diffArray(1,indexB6) = AngularDiff(centcoord4(indexB5,:),xyz4(indexB6,:));
                        end
                        diff1(1,indexB5) = min(diffArray);
                    end
                    diffsum = sum(diff1);
                    if diffsum<minDiff
                        minDiff = diffsum;
                        minOrient(1,1) = indexB7;
                        minOrient(1,2) = indexB8;
                        minOrient(1,3) = indexB9;
                    end
            end
        end
    end
    % Orientations in reverse order for proper intrinsic angles
    [phi1]=minOrient(1,3);
    [theta]=minOrient(1,2);
    [phi2]=minOrient(1,1);
    
end
end
