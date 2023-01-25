function [centcoord]=splitpeaks(X)
% To be called by 'peakfinder' function if statement. Where there are 
%peaks at high > 330 and low < 30 azimuth, these are likely the same peak. 
%'splitpeaks' calculates the mean coordinate and removes the duplicate.
    centcoord1=X;
    startx=(360+centcoord1(1,1));  
    centcoord1(1,1)=(startx+centcoord1(end,1))/2;
        if centcoord1(1,1)>360
                centcoord1(1,1)=centcoord1(1,1)-360;
        end
    centcoord1(1,2) = ((centcoord1(1,2)+centcoord1(end,2))/2);    
    centcoord1(length(centcoord1),:)=[];
    centcoord=centcoord1;
end