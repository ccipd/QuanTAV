function [T_xy,T_xz,T_yz] = compute_cartesian_orientations(vessel_points, ...
                                         seg_t, pixel_dimensions, options)
%% COMPUTE_CARTESIAN_ORIENTATIONS
% Generate three 2D projections of the tumor vessels by flattening tumor
% vessels along each cartesian dimension (X, Y, and Z) and calculate a 
% distribution of vessel orientations along each. Within a sliding window 
% the top five most prominent vessel orientations are computed via Hough
% transform and stored. Outputs the accumulated vessel orientations 
% derived from XY, XZ, and YZ projections
%
% INPUT ARGUMENTS
%   vessel_points - cell of all coordinates of voxels belonging to 
%       vasculature 
%   seg_t - binary volume of the tumor 
%   pixel_dimensions - three element vector of X,Y,Z pixel spacing of input
%       volume. We observed best results when resizing all volumes to 
%       isotropic 1 mm resolution prior to feature extraction, but
%       differing resolutions between scans can still be accounted for with 
%       this setting
%   options - optional struct including more advanced settings. Some of 
%       these settings were tuned on a per-modality basis in cross-
%       validation, while others were set to a single, fixed value for all 
%       experiments
%    --Tuned Parameters--
%       options.distance_range - range of distances (in mm) from tumor of 
%           vessels used to create projection images. Set from 0 to 11 mm 
%           for breast MRI and 0 to 7 mm for lung CT. Default: [0,11]
%       options.window_size - width (in pixels) of square window used to 
%           compute local vessel orientations within each projection image.
%           Set to 35 pixels for breast MRI and 20 pixels for lung CT. 
%           Default: 35 pixels 
%    --Other Parameters--
%       options.window_step - step size (in pixels) for moving the sliding 
%           window across projection images. Was set automatically to 1/3
%           options.window_size and this is the default behaviour, but can
%           be overridden manually. 

if ~exist('options','var')
    options = struct;
end
options = parse_options(options); 

% flattens vessels along each image axis
[BW_xy, BW_xz, BW_yz] = create_2D_cartesian_projections(vessel_points, seg_t, ...
                              options.distance_range, pixel_dimensions); 

% Compute vessel orientations across each projection
[T_xy,~] = process_hough_2D(BW_xy, options.window_size, options.window_step, 0);
[T_xz,~] = process_hough_2D(BW_xz, options.window_size, options.window_step, 0);
[T_yz,~] = process_hough_2D(BW_yz, options.window_size, options.window_step, 0);

end


% assign default values for any unspecified settings
function options = parse_options(options)

if isfield(options,'distance_range')
    if length(options.distance_range) == 1
        options.distance_range = [0,options.distance_range];
    end
else
    options.distance_range = [0,11];
end

if ~isfield(options,'window_size')
    options.window_size = 35;
end

if ~isfield(options,'window_step')
    options.window_step = round(options.window_size/3);
end

end
