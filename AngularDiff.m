function [Error] = AngularDiff(VecC,VecTP) 
%spherical coordinates of the topo data (vecc) and model output (vectp). 
% Must be in radians!   
% Returns the angular difference between a topography vector and orthogonal
% 5-vector set by Ralph's formula. To be called by bestfit1 and 
% vectormatching functions.

            [XC,YC,ZC] = sph2cart(VecC(1,1),VecC(1,2),VecC(1,3));
            [XTP,YTP,ZTP] = sph2cart(VecTP(1,1),VecTP(1,2),VecTP(1,3));
            XYZC = [XC,YC,ZC];
            XYZTP = [XTP,YTP,ZTP];
            Error = rad2deg(acos(dot(XYZC,XYZTP)));
            Error = real(Error);
end