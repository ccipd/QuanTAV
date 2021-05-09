function [BW_phitheta,BW_rphi,BW_rtheta] = create_2D_spherical_projections(S_mm, ...
                                            seg_t, pixel_dimensions, distance_range, projection_size)

% Convert skeleton to spherical coordinates
[theta, phi, r] = cartesian_to_spherical(S_mm, seg_t, pixel_dimensions); 

% rescale rotation from [-pi,pi] to [1,projection_size]
theta = 1+round((.5*projection_size*(theta+pi)./pi));
% rescale elevation from [-pi/2,pi/2] to [1,projection_size]
phi = 1+round((.5*projection_size*(phi+pi/2)./(pi/2)));
% rescale distance from distance_range to [1,projection_size]
r_width = distance_range(2)-distance_range(1);
r = 1+round((r - distance_range(1)).*projection_size./r_width);

pt_mm = cell2mat(S_mm(:)); % flatten all branches into single matrix
% restore pixel coordinates to compare with tumor binary mask
pt_idx = pt_mm ./ pixel_dimensions; 

% Initialize projection images
BW_phitheta = zeros(projection_size+1,projection_size+1);
BW_rphi = zeros(projection_size+1,projection_size+1);
BW_rtheta = zeros(projection_size+1,projection_size+1);

% Iterate through points and construct projections 
for ii =1:length(r)
    % Only include points within distance range and outside the tumor
    if r(ii)<=(projection_size) && r(ii)>0 && ... 
            seg_t(round(pt_idx(ii,1)),max([round(pt_idx(ii,2)),1]),round(pt_idx(ii,3))) == 0 
        BW_phitheta(phi(ii),theta(ii)) = 1;
        BW_rphi(r(ii),phi(ii)) = 1;
        BW_rtheta(r(ii),theta(ii)) = 1;
    end
end

% reorient projections for better continuity/interpretability 
BW_phitheta = ...
    BW_phitheta(:,[round(.75*projection_size)+1:projection_size+1,1:round(.75*projection_size)]);
BW_rphi = BW_rphi';
BW_rtheta = ...
    BW_rtheta(:,[round(.75*projection_size)+1:projection_size+1,1:round(.75*projection_size)]);
BW_rtheta = BW_rtheta';

