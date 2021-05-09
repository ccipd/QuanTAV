function [ T_all,R_all ] = process_hough_2D(BW_image, Wn, step, pad_dim)
%% PROCESS_HOUGH_2D
% Function to compute the distribution of vessel orientations within a 2D
% projection image. A sliding window of size Wn is slid across an image
% at a spacing of size step. For each window including vessel pixels, the
% theta and rho values associated with the top 5 peaks in the Hough
% tranform matrix are computed and accumulated. T_all is the set of all
% orientations (given by theta) across the projection image.
%
% Specify any dimensions that are continuous (e.g. rotation, where 0 and 
% 359 degrees are adjacent but on opposite ends of the projection image) 
% with pad_dim to wrap the image with a pad of size Wn at that axis to 
% prevent discontinuity
%%Prateek Feb15-2018

% Pad a continuous dimension (e.g. 
if pad_dim ~= 0
    rep_image = repmat(BW_image, double(pad_dim==1) + 1,  double(pad_dim==2) + 1); 
    new_image = rep_image(1:(size(BW_image,1)+(pad_dim==1)*(Wn-1)), ...
        1:(size(BW_image,2)+(pad_dim==2)*(Wn-1))); 
    BW_image = new_image; 
end


s=size(BW_image);
T_all=[]; R_all=[]; % initialize accumulator variables

% slide window over image and compute local vessel orientations 
for i=1:step:s(1)-Wn
    for j=1:step:s(2)-Wn
        BW_roi=BW_image(i:i+Wn-1,j:j+Wn-1);
        % only compute orientations if window contains vessel
        if sum(BW_roi(:))>0
            [H,~,~] = hough(BW_roi,'thetaresolution',.5);
            P = houghpeaks(H,5);
            if nnz(H)>4
                R_all=[R_all P(:,1)'];
                T_all=[T_all P(:,2)'];
            end
        end
    end
end

% % If image smaller than window, compute one set of peaks over full image
% if min(s) < Wn || isempty(T_all)
%     BW_roi = BW_image; 
%     [H,~,~] = hough(BW_roi);
%     P = houghpeaks(H,5);
%     R_all=[R_all P(:,1)'];
%     T_all=[T_all P(:,2)'];
% end

end