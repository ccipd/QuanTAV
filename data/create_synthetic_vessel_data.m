function [seg_t,seg_v] = create_synthetic_vessel_data(n_tumor_vessels, n_other_vessels, smoothness, numpoints, sz, tumor_radius, main_vessel_draw_radius,max_num_branches,branch_probability,branching_radius, vessel_width)

if length(sz) == 1
    sz = [sz,sz,sz];
end
center = round(sz/2);
if length(smoothness) == 1
    smoothness = [smoothness, smoothness];
end
% if length(numpoints) == 1
%     numpoints = [numpoints, numpoints];
% end
if length(vessel_width) == 1
    vessel_width = [vessel_width, vessel_width];
end
if main_vessel_draw_radius < tumor_radius 
    main_vessel_draw_radius = Inf;
end

random_smoothness = @() smoothness(1) + rand()*(smoothness(2)-smoothness(1));
random_numpoints = @() randi(numpoints, 1);
random_width = @() randi(vessel_width, 1);
euclidean_distance = @(x,y,z) sqrt(sum(([x(1),y(1),z(1)]-[x(2),y(2),z(2)]).^2,2));


%% generate synthetic tumor as sphere centered within volume
origin = zeros(sz); 
origin(center(1),center(2),center(3)) = 1; 
origin_dist = bwdist(origin);
seg_t = bwdist(origin) < tumor_radius;

% coordinates of all points within tumor, for main branch starting points
[x_t, y_t, z_t] = ind2sub(sz,find(seg_t)); 
n_t = length(x_t);

% coordinates of all points outside tumor, for main branch ending points
[y_p, x_p, z_p] = ind2sub(sz,find(~seg_t & origin_dist<main_vessel_draw_radius)); 
n_p = length(x_p);

n_vessels = n_tumor_vessels + n_other_vessels;

S = {}; 
for ii = 1:n_vessels
    N = random_numpoints();

    
    % choose random starting point in tumor and random end point outside it
    if ii <= n_tumor_vessels
        x_st = x_t; y_st = y_t; z_st = z_t;
        i_st = randi(n_t, 1);
    else
        x_st = x_p; y_st = y_p; z_st = z_p;
        i_st = randi(n_p, 1);
    end
    i_p = randi(n_p, 1);
    
    d_tp = euclidean_distance([x_st(i_st),x_p(i_p)],[y_st(i_st),y_p(i_p)],[z_st(i_st),z_p(i_p)]);
    N=round(numpoints*d_tp);
    X = make_random_1D_line([x_st(i_st),x_p(i_p)],N,random_smoothness(),sz(1));
    Y = make_random_1D_line([y_st(i_st),y_p(i_p)],N,random_smoothness(),sz(2));
    Z = make_random_1D_line([z_st(i_st),z_p(i_p)],N,random_smoothness(),sz(3));

    
    S{end+1} = [Y,X,Z]; % account for differences in plotting vs. indexing axis order

    if rand() <= branch_probability
        num_branches = randi(max_num_branches); 
        for jj = 1:num_branches
            % choose random point on current vessel as starting point
            i_v = randi(length(X)); 
     
            
            d_b = sqrt(sum(([x_p,y_p,z_p]-[X(i_v),Y(i_v),Z(i_v)]).^2,2));
            eligible_branch_ends = d_b < branching_radius; 
            x_b = x_p(eligible_branch_ends);
            y_b = y_p(eligible_branch_ends);
            z_b = z_p(eligible_branch_ends);
            i_b = randi(length(x_b)); 
                        
            x2 = make_random_1D_line([X(i_v),x_b(i_b)],N,random_smoothness(),sz(1));
            y2 = make_random_1D_line([Y(i_v),y_b(i_b)],N,random_smoothness(),sz(2));
            z2 = make_random_1D_line([Z(i_v),z_b(i_b)],N,random_smoothness(),sz(3));
            
            X = [X;x2]; Y = [Y;y2]; Z = [Z;z2];
            S{end+1} = [y2,x2,z2]; % account for differences in plotting vs. indexing axis order
        end
    end
end
seg_v = zeros(size(seg_t)); 

for branch = 1:length(S)
    pt = S{branch}; 
    seg_branch = zeros(size(seg_v)); 
    for ii = 1:length(pt)
        xyz = round(pt(ii,:)); 
        xyz = max([xyz; ones(1,3)]); 
        xyz = min([xyz; size(seg_t)]); 
        if seg_t(xyz(1),xyz(2),xyz(3))==0 && all(xyz) >= 1 && all(xyz <= sz)
            seg_branch(xyz(1),xyz(2),xyz(3))=1;
        end
    end
    % dilate vessel a random amount
    v_width = random_width();
    if v_width > 0
        seg_branch = imdilate(seg_branch,strel('sphere',v_width)); 
    end
    seg_v = seg_v | seg_branch; 
end