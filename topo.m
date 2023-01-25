function [phi1,theta,phi2]=topo(surfFOV)
    %Called by CrysChemMain
    %For a given surface area, returns the centroid coordinates 
    % corresponding to the normal etch vectors. 
    
    % Smoothing the data but ignoring the nans
    nanMask = isnan(surfFOV);                     %Create mask
    [r, c] = find(~nanMask);                %Index the nan locations
    [rNan, cNan] = find(nanMask);           %v.s.
    S = scatteredInterpolant(c, r, surfFOV(~nanMask), 'nearest');  %Interpolate the nans over surface
    interpVals = S(cNan, rNan);
    z1 = surfFOV;                           %To save original data
    z1(nanMask) = interpVals;
    zsmooth = imgaussfilt(z1,1);            %Smoothing the data
    zsmooth(nanMask) = nan;                 %Replace nans into data :)
    
    % Extract gradients and directions
    [Gx, Gy] = imgradientxy(zsmooth);
    [Gmag,Gdir] = imgradient(Gx, Gy);       %For gradient magnitude and direction
    mAngle = atand(Gmag);                   %The slope angle from horizontal (degrees)
    Gdir5 = round(Gdir);                    %Round for precision error               
    mAngle2 = round(mAngle);                %Round for precision error
    Gdir6 = reshape(Gdir5,[length(Gdir5).*length(Gdir5),1]);
    mAngle3 = reshape(mAngle2,[length(mAngle2).*length(mAngle2),1]); 
    GdirmAngle1 = [Gdir6,mAngle3];
    
    % Intensity for cartesian plots
    [~,~,ic1] = unique(GdirmAngle1, 'rows', 'stable');   %Unique values retaining order 
    h1 = accumarray(ic1, 1);                              %Count occurrences
    maph1 = h1(ic1);                                       %Map to 'ic' values
    result1 = [GdirmAngle1, maph1];
    resultun1 = rmmissing(unique(result1,"rows"));        %Remove nan data and duplicates
    
    % Gridding for cartesian plots
    [xq1,yq1] = meshgrid(-180:1:179, 0:1:89);           %Query points for interp: Az & El (cartesian)
    vq1 = griddata(resultun1(:,1),resultun1(:,2),resultun1(:,3),xq1,yq1);     %Interp
    
    %peak find struggles with low elevation peaks associated with Al-(100)
    vqnan = vq1;
    vqnan(isnan(vqnan))=0;
    q=sum(vqnan(:,1:end),2);
    v=cumtrapz(q);
    centcoord=peakfinder(vqnan);

    if prctile(v,20)/max(v) > 0.2
        [phi1,theta,phi2]=vectormatching(centcoord);
    else
        [phi1,theta,phi2]=bestfit1(centcoord);
    end
end