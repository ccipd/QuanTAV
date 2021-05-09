function [BW_xy,BW_xz,BW_yz] = create_2D_cartesian_projections(vessel_points, ...
                                   seg_t, distance_range, pixel_dimensions)
% Creates 2D cartesian projection images from coordinates of voxels of
% vessel volume. Flattens vasculature along one (X,Y,Z) dimension at a time 
% to obtain images of vessel position along two axes (e.g. XY, XZ, YZ).
%
% INPUT ARGUMENTS
%   vessel_points - cell of all coordinates of voxels belonging to 
%       vasculature 
%   seg_t - binary volume of the tumor 
%   distance_range - vector of two values indicating minimum and maximum 
%       distance (in mm) from tumor surface for a point in the vessel 
%       volume to be included in projection images.
%   pixel_dimensions - three element vector of X,Y,Z pixel spacing of input
%       volume. 

% Convert vessel coordinates to physical spacing
vessel_points_mm = convert_skeletons_to_mm_and_interpolate(vessel_points, ...
                                                    pixel_dimensions, 0); 

% Compute distance from tumor surface
[~, ~, d] = cartesian_to_spherical(vessel_points_mm, seg_t, pixel_dimensions);

% Reconstruct 3D vessel volume, including only points within distance_range
% and outside the tumor volume
BW_xyz = zeros(size(seg_t));
pt = cell2mat(vessel_points(:));
for ii = 1:length(pt)
    xyz = round(pt(ii,:)); 
    xyz = max([xyz; ones(1,3)]); 
    xyz = min([xyz; size(seg_t)]); 
    d_pt = d(ii); 
    if seg_t(xyz(1),xyz(2),xyz(3))==0 && ...
            (d_pt > distance_range(1) && d_pt <=distance_range(2))
        BW_xyz(xyz(1),xyz(2),xyz(3))=1;
    end
end

% Take projections across each cartesian dimension
BW_xy = max(BW_xyz,[],3);
BW_xz = squeeze(max(BW_xyz,[],1));
BW_yz = squeeze(max(BW_xyz,[],2));
    
end

