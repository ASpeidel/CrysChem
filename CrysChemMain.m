%% CrysChemMain
% CrysChemMain takes sequentially named topography datasets (.xyz) and
% outputs variable globalStore, which is a 3D array of Euler angles Phi1, 
% theta, Phi2 for each xy position for each measurement space.

% The topography data here is made up of separate patterns acquired by the 
% inteferometer, which are acquired in rows of 20 discrete patterns.

%Oversampling is used to enhance mapping resolution, using a moving
%sampling window. This leads to undersampling of the bottom right, so a
%jigsawing approach is used to resample undersampled locations.

myFolder = 'C:\Topography_file_location'; % Specify the folder 
filePattern = fullfile(myFolder, '*.xyz'); 
theFiles = dir(filePattern);
globalStore=[];
mapTstore=[];

Cs=1;
Ce=20;
while Ce<=length(theFiles)
    for C=Cs:Ce
        fullFileName = fullfile(myFolder, theFiles(C).name) 
        [X, Y, Z] = importfiletest(fullFileName);
        z = xyz2grid(X,Y,Z);
        sz=size(z);
        boxoffsetR=round((sz(1,1)-(637*4))/2)+1;
        boxoffsetC=round((sz(1,2)-(637*4))/2)+1;
        zcrop=z(boxoffsetR:boxoffsetR+(4*637),boxoffsetC:boxoffsetC+(4*637));
        
        MapT=zeros(6,6,3);
        for index1=1:6 % 637 pixel 110 um step size
            for index2=1:6
                surfFOV=imcrop(zcrop,[(index1*319-318) (index2*319-318) 637 637]);
                [xi,yi]=meshgrid(1:length(surfFOV),1:length(surfFOV));
                [xo,yo,zo] = prepareSurfaceData(xi,yi,surfFOV);
                sf = fit([xo, yo],zo,'poly33'); % A 3rd order polynomial to remove surface form
                zo=feval(sf,[xo,yo]);
                zo=flipud(xyz2grid(xo,yo,zo));
                z3=surfFOV-zo;
                [phi1,theta,phi2]=topo(z3);
                MapT(index2,index1,1)=phi1;
                MapT(index2,index1,2)=theta;
                MapT(index2,index1,3)=phi2;
            end
        end
        
        % Jigsawing to improve resolution
        %Extract final column, row, and corner FOVs and manipulate so that
        %undersampled topography is at the top left
        MapCol=zeros(2,6,3);
        OSCol=fliplr(zcrop(:,end-(2*637):end)).';
        OSRow=rot90(zcrop(end-(2*637):end,:),2);
        OSCnr=rot90(zcrop(end-(2*637):end,1:(2*637)),-1);
        
        for index1=1:6
            for index2=1:2
                surfFOV=imcrop(OSCol,[(index1*319-318) (index2*319-318) 637 637]);
                surfFOV=flipud(surfFOV).';
                [xi,yi]=meshgrid(1:length(surfFOV),1:length(surfFOV));
                [xo,yo,zo] = prepareSurfaceData(xi,yi,surfFOV);
                sf = fit([xo, yo],zo,'poly33'); % A 3rd order polynomial to remove surface form
                zo=feval(sf,[xo,yo]);
                zo=flipud(xyz2grid(xo,yo,zo));
                z3=surfFOV-zo;
                [phi1,theta,phi2]=topo(z3);
                MapCol(index2,index1,1)=phi1;
                MapCol(index2,index1,2)=theta;
                MapCol(index2,index1,3)=phi2;
            end
        end
        MapCol1=fliplr(pagetranspose(MapCol));
                
        MapRow=zeros(2,6,3);
        for index1=1:6
            for index2=1:2
                surfFOV=imcrop(OSRow,[(index1*319-318) (index2*319-318) 637 637]);
                surfFOV=rot90(surfFOV,2);
                [xi,yi]=meshgrid(1:length(surfFOV),1:length(surfFOV));
                [xo,yo,zo] = prepareSurfaceData(xi,yi,surfFOV);
                sf = fit([xo, yo],zo,'poly33'); % A 3rd order polynomial to remove surface form
                zo=feval(sf,[xo,yo]);
                zo=flipud(xyz2grid(xo,yo,zo));
                z3=surfFOV-zo;
                [phi1,theta,phi2]=topo(z3);
                MapRow(index2,index1,1)=phi1;
                MapRow(index2,index1,2)=theta;
                MapRow(index2,index1,3)=phi2;
            end
        end
        MapRow1=rot90(MapRow,2); %Reverse the original manipulation
        MapCnr=zeros(2,2,3);
        
        for index1=1:2
            for index2=1:2
                surfFOV=imcrop(OSCnr,[(index1*319-318) (index2*319-318) 637 637]);
                surfFOV=rot90(surfFOV);
                [xi,yi]=meshgrid(1:length(surfFOV),1:length(surfFOV));
                [xo,yo,zo] = prepareSurfaceData(xi,yi,surfFOV);
                sf = fit([xo, yo],zo,'poly33'); % A 3rd order polynomial to remove surface form
                zo=feval(sf,[xo,yo]);
                zo=flipud(xyz2grid(xo,yo,zo));
                z3=surfFOV-zo;
                [phi1,theta,phi2]=topo(z3);
                MapCnr(index2,index1,1)=phi1;
                MapCnr(index2,index1,2)=theta;
                MapCnr(index2,index1,3)=phi2;
            end
        end
        MapCnr1=rot90(MapCnr);
        MapBotRow=[MapCnr1, MapRow1];
        MapT1=[MapT,MapCol1];
        MapT2=[MapT1;MapBotRow];
        mapTstore=[mapTstore MapT2];
                
    end
Cs=Cs+20;
Ce=Ce+20;
globalStore=[globalStore; mapTstore];
mapTstore=[];
end
toc
