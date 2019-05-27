% items:    Data matrix. Row: Item. Columns: Weight, Profit
% mode_ : FeatureIDs that will be calculated
% varargin: Transformation data (M, W) and kp info
function feature_ = getKPFeatures(items, mode_, varargin)
nbFeatures = size(mode_);
maxValue = max(items,[],1);
feature_ = zeros(nbFeatures);

isTransformed = false; M = -1; W = -1;
if ~isempty(varargin), isTransformed = true; M = varargin{1}; W = varargin{2}; end

meanData = mean(items,1) ./ maxValue;
medianData = median(items,1) ./ maxValue;
stdData = std(items,0,1) ./ maxValue;
corrData = corr(items(:,1), items(:,2));

% Original features
feature_(mode_ == 1) = meanData(1); % Normalized mean weight
feature_(mode_ == 2) = medianData(1);      % Normalized median weight
feature_(mode_ == 3) = stdData(1);     % Normalized standard deviation (weight)
feature_(mode_ == 4) = meanData(2);      % Normalized mean profit
feature_(mode_ == 5) = medianData(2);      % Normalized median profit
feature_(mode_ == 6) = stdData(2); % Normalized standard deviation (profit)
feature_(mode_ == 7) = corrData /2 + 0.5; % Correlation (weight - profit)
% pw = items(:,2)./items(:,1); feature_(mode_ == 8) = mean(pw)/max(pw); % Correlation (weight - profit)
% kp = varargin{3}; feature_(mode_ == 100) = kp.freeCapacity / kp.maxWeight; % Remaining capacity (normalized)
%                     packedItems = kp.items; if isempty(packedItems), feature_(mode_ == 101) = 1; else, feature_(mode_ == 101) = mean(packedItems(:,1)) / max(packedItems(:,1)); end % Normalized mean weight of packed items
%                                              if isempty(packedItems), feature_(mode_ == 102) = 1; else, feature_(mode_ == 102) = mean(packedItems(:,2)) / max(packedItems(:,2)); end % Normalized mean profit of packed items

if isTransformed
    % Linear transformed features
    validID = mode_ == 8; feature_(validID) = linearTransform( meanData(1), M(validID), W(validID) ); % Normalized mean weight
    validID = mode_ == 9; feature_(validID) = linearTransform( medianData(1), M(validID), W(validID) );      % Normalized median weight
    validID = mode_ == 10; feature_(validID) = linearTransform( stdData(1), M(validID), W(validID) );     % Normalized standard deviation (weight)
    validID = mode_ == 11; feature_(validID) = linearTransform( meanData(2), M(validID), W(validID) );      % Normalized mean profit
    validID = mode_ == 12; feature_(validID) = linearTransform( medianData(2), M(validID), W(validID) );      % Normalized median profit
    validID = mode_ == 13; feature_(validID) = linearTransform( stdData(2), M(validID), W(validID) ); % Normalized standard deviation (profit)
    validID = mode_ == 14; feature_(validID) = linearTransform( corrData /2 + 0.5, M(validID), W(validID) ); % Correlation (weight - profit)
    
    % S-shaped transformed features
    validID = mode_ == 15; feature_(validID) = exponentialTransform( meanData(1), M(validID), W(validID) ); % Normalized mean weight
    validID = mode_ == 16; feature_(validID) = exponentialTransform( medianData(1), M(validID), W(validID) );      % Normalized median weight
    validID = mode_ == 17; feature_(validID) = exponentialTransform( stdData(1), M(validID), W(validID) );     % Normalized standard deviation (weight)
    validID = mode_ == 18; feature_(validID) = exponentialTransform( meanData(2), M(validID), W(validID) );      % Normalized mean profit
    validID = mode_ == 19; feature_(validID) = exponentialTransform( medianData(2), M(validID), W(validID) );      % Normalized median profit
    validID = mode_ == 20; feature_(validID) = exponentialTransform( stdData(2), M(validID), W(validID) ); % Normalized standard deviation (profit)
    validID = mode_ == 21; feature_(validID) = exponentialTransform( corrData /2 + 0.5, M(validID), W(validID) ); % Correlation (weight - profit)
end






