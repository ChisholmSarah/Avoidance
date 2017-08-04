%% This function divides the time series into sub-intervals from midday to
% midday and then permutes these sub-intervals perm times. For each of the
% permutations the number of times the diad is in the Distances is counted.

function [ShiftIdx, pLess, pMore, CountLess, CountSame, CountMore, DataMatrixPerm, PermCount, DistPerm] = PermutationFun(LinIntData, LinIntArray, perm, Dist, DistCountLinInt)
    
    n               = size(LinIntArray);
    PermData        = repmat(LinIntArray, [1, 1, 1, perm]);
    Permutations    = NaN(1,perm);
    PermCount       = NaN(size(Dist,2)+1,perm);
    DataMatrixPerm  = [];
    DistPerm        = NaN(size(LinIntArray,1),perm);

    %% Note down the beginning index of each Block for LinIntArray
    ShiftIdx = size(LinIntData(LinIntData(:,5,1,1)~=0,:,1,1), 1);
    for i = 2:size(LinIntData,4) % dim 3 of LinIntData is per individual and dim 4 per block; start at 2, because we know the first one is 1
        ShiftIdx(i,1) = ShiftIdx(i-1,1)+size(LinIntData(LinIntData(:,5,1,i)~=0,:,1,i), 1); % row numbers of the last measurements of each block
    end

    %% Create perm permutations, and add a column to the first
    % matrix (only one of the individuals needs to be permuted)
    % which has a number per block. The matrices are then sorted
    % first according to that column, then according to the
    % cumulative time column (the 6th column).

    if ~isempty(ShiftIdx)
        for p=1:perm
            s                                  = RandStream('mt19937ar','Seed',p);
            Permutations(1:size(ShiftIdx,1),p) = randperm(s, size(ShiftIdx,1));
            % Add a column to LinIntArray with the permuted index per
            % block
            PermData(1:ShiftIdx(1),n(2)+1,1,p) = ones(ShiftIdx(1),1)*Permutations(1,p);
            for j=2:length(ShiftIdx)
                PermData((ShiftIdx(j-1)+1):ShiftIdx(j),n(2)+1,1,p) = ones(ShiftIdx(j)-ShiftIdx(j-1),1)*Permutations(j,p);
            end
            x                           = PermData(:,:,1,p);
            x                           = x(x(:,5)~=0,:);
            PermData(1:size(x,1),:,1,p) = sortrows(x, [n(2)+1, n(2)]);

            %% Calculate the distance between the permuted and the non-permuted individuals
            x                         = PermData(:,:,:,p);
            DataMatrixPerm            = x(x(:,5,1)~=0,:,:);
            [x1, x2]                  = DistanceFun(DataMatrixPerm, Dist);
            DistPerm(1:size(x1,1),p)  = x1;
            PermCount(1:size(x2,1),p) = x2;
        end
    end
    
    if ~isempty(PermCount)
        LinIntCount = repmat(DistCountLinInt, 1, perm);
        Less        = PermCount<LinIntCount;
        Same        = PermCount==LinIntCount;
        More        = PermCount>LinIntCount;
        CountLess   = sum(Less,2); % Number of times the Permuted data is less often within that distance of each other than the observed
        CountSame   = sum(Same,2);
        CountMore   = sum(More,2);

        %% If pLess<0.05/(size(Dist,2)*2) for a particular distance, then that
        % dyad is significantly less often in that particular distance
        % than expected by chance. Similarly if
        % pMore<0.05/(size(Dist,2)*2) then the dyad is significantly
        % more often in that distance than expected.

        pLess       = (CountLess+CountSame)/perm;
        pMore       = (CountMore+CountSame)/perm;
    else
        CountLess   = [];
        CountSame   = [];
        CountMore   = [];
        pLess       = [];
        pMore       = [];
    end
    
end
