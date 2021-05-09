function stats = quickstats(feats)
stats = [];
if strcmp(class(feats), 'cell')
    for i = 1:length(feats)
        x = feats{i};
        if size(x,1) == 1
            x = repmat(x,2,1);
        end
        
        stats = [stats,mean(x,1),median(x,1),std(x,1),skewness(x,1),kurtosis(x,1)];
    end
else
    stats = [mean(feats),median(feats),std(feats),skewness(feats),kurtosis(feats)];
end
end