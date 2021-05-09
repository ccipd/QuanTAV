function [T_phitheta, T_rtheta, T_rphi] = compute_spherical_orientations(S, ...
                                            seg_t, pixel_dimensions, options)
%% COMPUTE_SPHERICAL_ORIENTATIONS
% Generate three 2D projections of the tumor vessels by flattening tumor
% vessels along each spherical dimension (distance from tumor surface, 
% rotation and elevation from the tumor centroid) and calculate a 
% distribution of vessel orientations along each. Within a sliding window 
% the top five most prominent vessel orientations are computed via Hough
% transform and stored. Outputs the accumulated vessel orientations 
% derived from elevation-rotation, distance-rotation, and distance-elevation 
% projections
%
% INPUT ARGUMENTS
%   S - skeleton of vessel network divided into branches. obtained from 
%       skeleton.m as:
%           S = skeleton(seg_v | seg_t); 
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
%       options.spherical_projection_size - width (in pixels) of square
%           spherical projection images generated from vessels. A spherical
%           dimension's range (distance_range for distance, 0-360 degrees
%           for rotation, -90-90 degrees for elevation) will be mapped
%           equally across this many pixels for each side of the projection
%           image. Fixed at default value of 400 pixels. 
%       options.n_interpolation - the number of interpolation points
%           inserted between adjacent points of vasculature to ensure 
%           continuity of spherical vessel projections. Set to fixed value
%           of 50 points, but does not have much bearing on results
%           provided set high enough for smooth vessel projections. 

if ~exist('options','var')
    options = struct;
end
options = parse_options(options); 

% Convert skeleton coordinates to physical distance and interpolate b/w
% points for continuous spherical projections
S_mm = convert_skeletons_to_mm_and_interpolate(S, pixel_dimensions, ...
                                                    options.n_interpolate); 

% Converts skeleton points to spherical coordinate system and projects 
% vessels along each spherical dimension
[BW_phitheta, BW_rphi, BW_rtheta] = create_2D_spherical_projections(S_mm, ...
                      seg_t, pixel_dimensions, options.distance_range, ...
                      options.spherical_projection_size); 

% Compute vessel orientations across each projection
[T_phitheta,~] = process_hough_2D(BW_phitheta,options.window_size,options.window_step,2);
[T_rtheta,~] = process_hough_2D(BW_rtheta,options.window_size,options.window_step,1);
[T_rphi,~] = process_hough_2D(BW_rphi,options.window_size,options.window_step,0);
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

if ~isfield(options,'spherical_projection_size')
    options.spherical_projection_size = 400;
end

if ~isfield(options,'n_interpolate')
    options.n_interpolate = 50;
end

end