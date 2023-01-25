function [BW1] = Imageproc(X)
% image area opening and closing to allow coordinate identification. 
% To be called by peakfinder function
BW1=X;
BW1=bwareaopen(BW1,10); 
se = strel('disk',10);
BW1=imclose(BW1,se);
end