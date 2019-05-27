% bins:   Structure with the following fields:
%             nbBins        = Current number of bins
%             bin           = Cell array with the information (items) of each bin
%             capacity      = Vector with maximum capacity of each bin
%             freeCapacity  = Vector with remaining capacity of each bin
%             defaultCapacity = Default value for all bins (assuming a single capacity for all)
% items: Remaining elements of the instance
% mode_ : FeatureIDs that will be calculated
function feature_ = getBinPackingFeatures(bins, items, mode_, varargin)

isPCA = false;

nbFeatures = size(mode_);
maxValue = ones(nbFeatures);
feature_ = zeros(nbFeatures);
isTransformed = false;

if ~isempty(varargin)
    if length(varargin) == 2
        isTransformed = true; M = varargin{1}; W = varargin{2};
    elseif length(varargin) == 3
        isPCA = true;
        mu = varargin{1}; sigma = varargin{2}; U = varargin{3};
    end
end


nbOpenBins = sum( bins.freeCapacity ~= 0 );
openRatio = nbOpenBins / bins.nbBins;
fullRatio = 1 - bins.freeCapacity ./ bins.capacity; % Fullness == 1 means perfect fit
avgFullness = 1/bins.nbBins * sum(fullRatio);
avgFullnessGen = 1/bins.nbBins * sum(fullRatio.^2);
avgPerfection = 1/bins.nbBins * sum(fullRatio(fullRatio==1));
avgLength = mean(items);
stdLength = std(items);
bigLocations = items > (bins.defaultCapacity/2); 
    bigPieces = sum(bigLocations) / length(items);
    bigItems = items(bigLocations);
        bigPieceAvg = mean(bigItems);
        bigPieceStd = std(bigItems);
smallLocations = items < (bins.defaultCapacity/4);
    smallPieces = sum(smallLocations) / length(items);
    smallItems = items(smallLocations);
        smallPieceAvg = mean(smallItems);
        smallPieceStd = std(smallItems);
notBigPieces = items(~bigLocations);
    notBigAvg = mean(notBigPieces);
    notBigStd = std(notBigPieces);

% Original features
feature_(mode_ == 1) = openRatio; % Ratio of open bins (original)
feature_(mode_ == 2) = avgFullness; % Average fullness
feature_(mode_ == 3) = avgFullnessGen; % Average generic fullness
feature_(mode_ == 4) = avgPerfection; % Average Perfection
feature_(mode_ == 5) = avgLength; % Average Length of remaining pieces
feature_(mode_ == 6) = stdLength; % Std dev of Length of remaining pieces
feature_(mode_ == 7) = bigPieces; % Ratio of big pieces
feature_(mode_ == 8) = bigPieceAvg; % Ratio of small pieces
feature_(mode_ == 9) = bigPieceStd; % Ratio of small pieces
feature_(mode_ == 10) = smallPieces; % Ratio of small pieces
feature_(mode_ == 11) = smallPieceAvg; % Ratio of small pieces
feature_(mode_ == 12) = smallPieceStd; % Ratio of small pieces
feature_(mode_ == 13) = notBigAvg; % Ratio of small pieces
feature_(mode_ == 14) = notBigStd; % Ratio of small pieces

totalFeatures = 14;

if isTransformed
    % Linear transformed features
    validID = mode_ == totalFeatures + 1; feature_(validID) = linearTransform( openRatio, M(validID), W(validID) ); % Ratio of open bins (transformed)
    validID = mode_ == totalFeatures + 2; feature_(validID) = linearTransform( avgFullness, M(validID), W(validID) ); % Average fullness (transformed)
    validID = mode_ == totalFeatures + 3; feature_(validID) = linearTransform( avgFullnessGen, M(validID), W(validID) ); % Average Generic fullness (transformed)
    validID = mode_ == totalFeatures + 4; feature_(validID) = linearTransform( avgPerfection, M(validID), W(validID) ); % Average Perfection (transformed)
    validID = mode_ == totalFeatures + 5; feature_(validID) = linearTransform( avgLength, M(validID), W(validID) ); % Average Perfection (transformed)
    validID = mode_ == totalFeatures + 6; feature_(validID) = linearTransform( stdLength, M(validID), W(validID) ); % Average Perfection (transformed)
    validID = mode_ == totalFeatures + 7; feature_(validID) = linearTransform( bigPieces, M(validID), W(validID) ); % Average Perfection (transformed)
    
    % S-shaped transformed features
    validID = mode_ == 2*totalFeatures + 1; feature_(validID) = exponentialTransform( openRatio, M(validID), W(validID) ); % Normalized mean weight
    validID = mode_ == 2*totalFeatures + 2; feature_(validID) = exponentialTransform( avgFullness, M(validID), W(validID) ); % Average fullness (transformed)
    validID = mode_ == 2*totalFeatures + 3; feature_(validID) = exponentialTransform( avgFullnessGen, M(validID), W(validID) ); % Average fullness (transformed)
    validID = mode_ == 2*totalFeatures + 4; feature_(validID) = exponentialTransform( avgPerfection, M(validID), W(validID) ); % Average fullness (transformed)
    validID = mode_ == 2*totalFeatures + 5; feature_(validID) = exponentialTransform( avgLength, M(validID), W(validID) ); % Average fullness (transformed)
    validID = mode_ == 2*totalFeatures + 6; feature_(validID) = exponentialTransform( stdLength, M(validID), W(validID) ); % Average fullness (transformed)
    validID = mode_ == 2*totalFeatures + 7; feature_(validID) = exponentialTransform( bigPieces, M(validID), W(validID) ); % Average fullness (transformed)
    
end

if isPCA
    feature_ = applyPCA(feature_, mu, sigma, U);
end

end