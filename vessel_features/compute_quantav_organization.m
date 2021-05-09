function organization_ftrs = compute_quantav_organization(S, seg_v, seg_t, pixel_dimensions, options)
%% COMPUTE_QUANTAV_ORGANIZATION
% Compute set of 30 quantitative tumor-associated vasculature (QuanTAV) 
% features describing the organization of vessels surrounding a tumor. 
% Six 2D projections of the tumor vessels are created along each cartesian
% (X, Y, and Z) and spherical (rotation and elevation from the tumor
% surface, distance from the tumor surface) dimension. A sliding window 
% is used to compute and  store local vessel orientations across each 
% vessel projection. Each set of five features consists of first order 
% statistics (mean, median, standard deviation, skewness, and kurtosis) for
% a projection image in the order XY, XZ, YZ, rotation-elevation, distance-
% rotation, and distance-elevation. 
%
% INPUT ARGUMENTS
%   S - skeleton of vessel network divided into branches. obtained from 
%       skeleton.m as:
%           S = skeleton(seg_v | seg_t); 
%   seg_v - binary volume of the tumor-adjacent vasculature
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
%          --Tuned Parameters--
%       options.distance_range - range of distances (in mm) from tumor of 
%           vessels used to create projection images. Set from 0 to 11 mm 
%           for breast MRI and 0 to 7 mm for lung CT. Default: [0,11]
%       options.window_size - width (in pixels) of square window used to 
%           compute local vessel orientations within each projection image.
%           Set to 35 pixels for breast MRI and 20 pixels for lung CT. 
%           Default: 35 pixels 
%          --Other Parameters--
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

% Get list of all points within tumor vasculature
c = find(seg_v);
[x, y, z] = ind2sub(size(seg_v),c);
vessel_points = {}; 
vessel_points{1} = [x,y,z];
% Create cartesians projections and compute orientations from all points in 
% vessel volume
[Txy, Txz, Tyz] = compute_cartesian_orientations(vessel_points, seg_t, ...
                                pixel_dimensions, options);

% Create spherical projections and compute orientations from vessel 
% centerline skeleton 
[Ttp, Trt, Trp] = compute_spherical_orientations(S, seg_t, ...
                                pixel_dimensions, options);
                            
cartesian_ftrs = [quickstats({Txy'}),quickstats({Txz'}),quickstats({Tyz'})];

                            
% Shift distance-angle projections to remove discontinuity for vessels 
% approaching horizontal
Trt(Trt>180) = Trt(Trt>180) - 360;
Trp(Trp>180) = Trp(Trp>180) - 360;

spherical_ftrs = [quickstats({Ttp'}),quickstats({Trt'}),quickstats({Trp'})];

organization_ftrs = [cartesian_ftrs, spherical_ftrs];

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