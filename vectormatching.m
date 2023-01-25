function [phi1,theta,phi2]=vectormatching(centcoord)
% To be called by 'topo'. Indexes the lattice around the low angle (high
% intensity) peak associated with crystals vicinal to (100)

if min(size(centcoord))<2 % Topography not indexable when only 1 centcoord peak is returned by peakfinder function
    [phi1]=NaN;
    [theta]=NaN;
    [phi2]=NaN;
else
xyz1 = [1 0 0; 0 1 0; 0 0 1;-1 0 0; 0 -1 0]; % Points to rotate
centcoord1=centcoord;
centcoord1(:,2)=90-centcoord(:,2);
centcoord1=deg2rad(centcoord1);
len=length(centcoord);
centcoord3=ones(len,1);
centcoord1=[centcoord1,centcoord3];

phi2diff=zeros(360,1);                  %The sum of the minimum separations
[~,I]=min(centcoord(:,2));             %Find the low angle coordinate
TF4 = zeros(360,1);                     %To find local minima
len=size(centcoord)-1;
Errors3d = zeros(360,4,len(1,1));       %To pre-allocate arrays for comparison
Errorsmin = zeros(360,1,len(1,1));

% Extract Euler angles phi1 and theta from topography data from xyz.m output (centcoord)     
R1 = (centcoord(I,1)); %Rotations about these angles match a model vector with the topography data output
R2 = (centcoord(I,2));

    for R3 = 1:360
        Rz1 = rotz(R3); %matrix rotations in reverse order for intrinsic Euler angles
        xyz2 = xyz1*Rz1;
        Ry = roty(-R2);
        xyz3 = xyz2*Ry;
        Rz2 = rotz(-R1);
        xyz4 = xyz3*Rz2;
        
        [az,el,r]=cart2sph(xyz4(:,1),xyz4(:,2),xyz4(:,3));
        xyz4 = horzcat(az,el,r);
        xyz4(:,1) = rad2deg(xyz4(:,1));
        centcoord4 = centcoord1;
    
        for index1 = 1:length(xyz4)
            if xyz4(index1,1)<0 
               xyz4(index1,1) = xyz4(index1,1)+360;  %Correcting for circular coordinate convention
            end
        end
        xyz4(:,1) = deg2rad(xyz4(:,1));
        %Routine to remove the matched vectors. Find and remove both instances
        UnqArray = [centcoord4;xyz4]; 
        UnqArray = round(UnqArray,4);   %Round for precision errors otherwise 'unique' function is temperamental
        [~,~,ic] = unique(UnqArray,'rows','stable');  % Find unique and return indices 
        flag=0;
        for i1 = 1:length(centcoord4)
               for i2 = length(centcoord4)+1:length(ic)
                    if ic(i1) == ic(i2)
                        flag=1;
                        break
                    end
               end
               if flag==1
                   break
               end
        end
        UnqArray(i1,:) = [];        %Delete first instance
        UnqArray(i2-1,:) = [];      %-1 as first instance is deleted in row above
        len=size(centcoord4);
        centcoord4 = UnqArray(1:len(1,1)-1,:);
        xyz5 = UnqArray(length(centcoord):length(UnqArray),:);
        len=size(centcoord4);       %Rewrite as length required for next loop
    
        for index3=1:len(1,1)
            for index5 = 1:length(xyz5)
                Errors3d(R3,index5,index3) = AngularDiff(centcoord4(index3,:),xyz5(index5,:));
            end
        end
        
        for index6 = 1:len(1,1)
            for index4 = 1:360
                Errorsmin(index4,1,index6) = min(Errors3d(index4,:,index6));
                phi2diff(index4,1) = sum(Errorsmin(index4,1,1:end));
            end
        end
    end

    %Find the local minima 
    TF4(:,1) = islocalmin(phi2diff(:,1),'MinProminence',1); %This cannot find minima crossing the azimuth origin 360 - 0
    if sum(TF4(:,1))>4 %To correct if this is the case, cut and paste last 45 deg to start and refind minima
       phisplit = phi2diff(:,1);
       phisplitb = phisplit(316:end,1);
       phisplit = phisplit(1:315,1);
       phisplit = [phisplitb;phisplit];
       TF4b = islocalmin(phisplit,'MinProminence',1);
       TF4split = TF4b(1:45,1);
       TF4b = TF4b(46:end,1);
       TF4b = [TF4b;TF4split];
       TF4(:,1) = TF4b;
    end
    
    if sum(TF4(:,1))<4 %To correct if this is the case, cut and paste last 45 deg to start and refind minima
       phisplit = phi2diff(:,1);
       phisplitb = phisplit(316:end,1);
       phisplit = phisplit(1:315,1);
       phisplit = [phisplitb;phisplit];
       TF4b = islocalmin(phisplit,'MinProminence',1);
       TF4split = TF4b(1:45,1);
       TF4b = TF4b(46:end,1);
       TF4b = [TF4b;TF4split];
       TF4(:,1) = TF4b;
    end

     %Index the local minima 
     TF5(:,1) = find(TF4(:,1));
     
     if length(TF5)<1
         [phi1,theta,phi2]=bestfit1(centcoord);
     else
     %In reverse order for itrinsic rotations
     [phi1]=TF5(1,1);
     [theta]=-round((centcoord(I,2)));
     [phi2]=-round((centcoord(I,1)));
     end
end