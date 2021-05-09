function [morphology_features,varargout] = compute_quantav_morphology(seg_v,seg_t,S)
%% COMPUTE_QUANTAV_ORGANIZATION
% Compute set of 61 quantitative tumor-associated vasculature (QuanTAV) 
% features describing the morphology of vessels surrounding a tumor through
% attributes such as torsion and curvature. 
%
% INPUT ARGUMENTS
%   seg_v - binary volume of the tumor-adjacent vasculature
%   seg_t - binary volume of the tumor
%   S - skeleton of vessel network divided into branches. obtained from 
%       skeleton.m as:
%           S = skeleton(seg_v | seg_t);  
%
% This set of QuanTAV Morphology features is expanded from the set of 35 
% proposed in "Quantitative vessel tortuosity: A potential CT imaging 
% biomarker for distinguishing lung granulomas from adenocarcinomas,"
% Alilou et al., Nature Sci Rep 2018. Original code, Mehdi Alilou 2018.
% Updated and expanded, Nathaniel Braman 2019.

% remove voxels within tumor from vessel segmentation
seg_v(seg_t==1) = 0;

% Initialize counters and accumulators 
sumL=0;
pt=1;
feeding_vessel_count = 0;
curvALL=zeros(1,10000); 

branches = {}; 

% Iterate through vascular skeleton one branch at a time  
for n_branch=1:length(S)
    L=S{n_branch}; % get branch

    %% filter points inside nodule
    pitem=0;
    L2=[0,0,0];
    enters_tumor = false;
    total_pts = 0; 
    branch_info = struct;
    
    pt_in_tumor = zeros(size(L,1),1); 
    
    for n_pt=1:size(L,1)
        % identify points inside the tumor
        if seg_t(round(L(n_pt,1)),round(L(n_pt,2)),round(L(n_pt,3))) ==1
            pt_in_tumor(n_pt) = 1; 
        end
    end
    % remove vessel points entering tumor
    L(pt_in_tumor==1,:) = []; 
    
    %% Remove adjacent points 
    % Sometimes the skeletonization code will place two points practically 
    % on top of each other, which causes curvature values to blow up
    % We compute the distance between adjacent points and delete the first
    % point in any pair separated by a distance < .001
    dists=sqrt((L(1:end-1,1)-L(2:end,1)).^2+(L(1:end-1,2)-L(2:end,2)).^2+(L(1:end-1,3)-L(2:end,3)).^2);
    pt_too_close = find(dists<.001) + 1;
    L(pt_too_close,:) = [];

    %% Compute branch-level vessel metrics
    if size(L,1)>10 % only analyze branches 11+ points long
        sumL=sumL+size(L,1);

        % compute euclidean distance from branch start and endpoints
        x=[L(1,1) L(size(L,1),1)];y=[L(1,2) L(size(L,1),2)];z=[L(1,3) L(size(L,1),3)];
        dist = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2+(z(1)-z(2))^2);
        
        % initialize length variable for measuring full branch length
        len = sqrt((L(1,1)-L(2,1))^2+(L(1,2)-L(2,2))^2+(L(1,3)-L(2,3))^2);

        curv=zeros(size(L,1),1);

        % iterate through branch and compute curvature for each set of 
        % 3 adjacent points
        for n_pt=2:size(L,1)-1
            p1=L(n_pt-1,:);p2=L(n_pt,:);p3=L(n_pt+1,:);
            
            % compute local curvature between p1, p2, p3
            [~,rad,~,~] = circlefit3d(p1,p2,p3);
            if isnan(rad)
                curv(n_pt,1) = 0;
            else
                curv(n_pt,1)=1/rad;
            end
            
            % add segment length to total branch length 
            len = len + sqrt((L(n_pt,1)-L(n_pt+1,1))^2 ...
                            +(L(n_pt,2)-L(n_pt+1,2))^2 ...
                            +(L(n_pt,3)-L(n_pt+1,3))^2);
        end
        curv(isnan(curv)) = [];
        
        % compute and store branch-level measures:
        %  - torsion of branch 
        %  - statistics of curvature measures across branch (mean, std, max, 
        %       skewness, and kurtosis)
        % + add branch curvature measures to set of all curvature measures
        if length(curv)>2 && length(unique(curv))>1
            % store branch-level metrics
            branch_torsion = (1-(dist/len));
            torsion(n_branch) = branch_torsion;
            curv_std(n_branch) = std(curv(2:size(curv,1)-1));
            curv_mean(n_branch) = mean(curv(2:size(curv,1)-1));
            curv_max(n_branch) = max(curv(2:size(curv,1)-1));
            curv_skew(n_branch) = skewness(curv(2:size(curv,1)-1));
            curv_kurt(n_branch) = kurtosis(curv(2:size(curv,1)-1));
            
            % accumulate curvature measures for full vessel network one 
            % branch at a time
            curv=curv';
            curvALL(1,pt:pt+(size(curv,2)-2)-1)=curv(1,2:size(curv,2)-1);
            pt=pt+size(curv,2)-2;

            % tally number of vessels that directly feed tumor
            enters_tumor = any(pt_in_tumor);
            if enters_tumor
                feeding_vessel_count = feeding_vessel_count +1;
            end
            
            branch_info.points = L; 
            branch_info.torsion = branch_torsion; 
            branch_info.curvature = curv;
            branch_info.enters_tumor = enters_tumor; 
            branches{end+1} = branch_info; 
        end

    end

end

% remove indices from branch-level measures corresponding to branches that 
% were skipped
curvALL=curvALL(1,1:min(pt,length(curvALL)));
torsion(torsion==0) = [];
curv_mean(curv_mean==0) = [];
curv_std(curv_std==0) = [];
curv_max(curv_max==0) = [];
curv_skew(curv_skew==0) = [];
curv_kurt(curv_kurt==0) = [];

% compute 10-bin histograms for global curvature and torsion measures
hiscurv=hist(curvALL);
histor=hist(torsion);

vessel_vol = sum(sum(sum(seg_v))); 
tumor_vol = sum(sum(sum(seg_t))); 

% % TODO: reorder
% morphology_features_ = [mean(torsion(torsion~=0)),...
%     std(torsion(torsion~=0)),max(torsion(torsion~=0)),...
%     mean(curv_std(curv_std~=0)),std(curv_std(curv_std~=0)),...
%     max(curv_std(curv_std~=0)),mean(curv_mean(curv_mean~=0)),...
%     std(curv_mean(curv_mean~=0)),max(curv_mean(curv_mean~=0)),...
%     mean(curv_max(curv_max~=0)),std(curv_max(curv_max~=0)),...
%     max(curv_max(curv_max~=0))  ,sumL,sum(sum(sum(seg_v)))/(size(seg_v,1)*size(seg_v,2)*size(seg_v,3)),...
%     sum(sum(sum(seg_v))),sum(sum(sum(seg_v)))/sum(sum(sum(seg_t))),hiscurv,histor, ...
%     feeding_vessel_count, feeding_vessel_count/length(torsion), ...
%     mean(curvALL),std(curvALL),max(curvALL),skewness(curvALL),...
%     kurtosis(curvALL),skewness(torsion),kurtosis(torsion),...
%     skewness(curv_mean),kurtosis(curv_mean),skewness(curv_std),...
%     kurtosis(curv_std),skewness(curv_max),kurtosis(curv_max),...
%     mean(curv_skew),std(curv_skew),max(curv_skew),skewness(curv_skew),...
%     kurtosis(curv_skew),mean(curv_kurt),std(curv_kurt),max(curv_kurt),...
%     skewness(curv_kurt),kurtosis(curv_kurt)];

% TODO: reorder
morphology_features = [vessel_stats(torsion), vessel_stats(curv_mean) ...
    vessel_stats(curv_std), vessel_stats(curv_max), vessel_stats(curv_skew), ...
    vessel_stats(curv_kurt), vessel_stats(curvALL),hiscurv,histor, ...
    vessel_vol, vessel_vol/(size(seg_v,1)*size(seg_v,2)*size(seg_v,3)), ... 
    vessel_vol/tumor_vol, sumL,feeding_vessel_count, ...
    feeding_vessel_count/length(torsion)];

if nargout > 1
    varargout{1} = branches;
end

end

function vessel_stats = vessel_stats(features)
vessel_stats = [mean(features),std(features),max(features),...
                        skewness(features),kurtosis(features)];
end
