function S_mm = convert_skeletons_to_mm_and_interpolate(S, px_dim, n_int)
%% CONVERT_SKELETONS_TO_MM_AND_INTERPOLATE
% Converts skeleton of pixel coordinates to physical spacing based on 
% acquisition resolution, px_dim, and interpolates between skeleton points 
% by adding n_int equally spaced coordinates between them. 
for ii = 1:length(S)
    S{ii} = S{ii}.*px_dim;
    S_mm{ii} = S{ii}(1,:); 
    for jj = 1:size(S{ii},1)-1
        pt1 = S{ii}(jj,:);
        pt2 = S{ii}(jj+1,:);
        int1 = linspace(pt1(1), pt2(1), n_int+2);
        int2 = linspace(pt1(2), pt2(2), n_int+2);
        int3 = linspace(pt1(3), pt2(3), n_int+2);
        S_mm{ii} = [S_mm{ii}; [int1(2:end)', int2(2:end)', int3(2:end)']];
    end
end

end
