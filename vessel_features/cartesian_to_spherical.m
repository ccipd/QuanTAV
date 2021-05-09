function [theta, phi, r] = cartesian_to_spherical(S_mm, seg_t, pixel_dimensions)
%% CARTESIAN_TO_SPHERICAL
% Convert skeleton coordinates (S_mm) with physical spacing in mm 
% (derived from pixel coordinates with spacing pixel_dimensions) to 
% spherical coordinate system corresponding to tumor segmentation (seg_t):
%   theta, rotation about the z axis through the tumor centroid 
%   phi, elevation from the XY plane at the tumor centroid 
%   r, the distance from the surface of the tumor volume 

% compute tumor centroid position, in pixels and mm
cent = regionprops3(seg_t,'Centroid');
[c] = cent.Centroid;
c_mm = c.*pixel_dimensions; 

pt_mm = cell2mat(S_mm(:)); % flatten all branches into single matrix
pt = pt_mm ./ pixel_dimensions;

% for each voxel in the vessel volume, get the corresponding position of 
% the nearest point on the tumor surface (x_s,y_s,z_s)
[~, idx_px] = bwdist(seg_t); 
[x_s, y_s, z_s] = ind2sub(size(seg_t), ...
        idx_px(sub2ind(size(seg_t),round(pt(:,1)), round(pt(:,2)),round(pt(:,3)))));

% compute r, the euclidean distance of each vessel point from tumor surface
nearest_surface_pt_mm = [x_s,y_s,z_s] .* pixel_dimensions; 
r = sqrt(sum((pt_mm  - nearest_surface_pt_mm).^2,2));

% compute spherical coordinates relative to tumor centroid
[theta,phi,~] = cart2sph(pt_mm(:,2)-c_mm(1), pt_mm(:,1)-c_mm(2), pt_mm(:,3)-c_mm(3));
end 