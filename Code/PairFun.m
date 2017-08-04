%% This function combines two files and divides the data into blocks
%  (Divided - one matrix per block). Blocks that start after half way
%  through the block, or end before half way through the block are
%  discarded. If both individuals are not tagged for Threshold hours,
%  discard.

% OUTPUT:   D       - array (kx6xn) with a matrix (:,:,i) per period
%           D1      - D matrix of the first input file (InputFile1)
%           D2      - D matrix of the second input file (InputFile2)
%           Lengths - a matrix (1,n) with the number of location fixes for
%               each of the n joint periods
%           Times   - an array (3x2xn) with the start time and date (first
%               row), the end time and date (second row) and the number of
%               days abd seconds (third row) for each of the n periods
%               (matrix per period)

function [D, D1, D2, DOld, Divided, Lengths, Times] = PairFun(D1, D2, Name1, Name2, Threshold, BlockStart, BlockLength)
    
    Divided = [];
    D       = [];
    Lengths = [];
    Times   = [];

    % if file 1 ends before file 2 starts or file 2 ends before file 1
    % starts, i.e. if they do not overlap, then stop and give warning
    % message.
    if (D1(end,size(D1,2))<=D2(1,size(D2,2)) || D2(end,size(D2,2))<=D1(1,size(D1,2)))
        warning('%s and %s are not tagged simultaneously', Name1, Name2);
    else
        
%% Find the overlapping observation periods and discard the rest
        DOld         = [D1; D2]; % join the two matrices into one
        DOld         = sortrows(DOld, 6); % sorting by the continuous time column
        Idx          = 1:(size(DOld,1)-1);
        DifferentInd = Idx(DOld(1:(end-1),5)~=DOld(2:end,5)); % vector of indices that indicate that two consecutive taggs are from different individuals 
        i            = DifferentInd(1,1); % i records the index of the first observation in each block
        
%% Divide the joint data into blocks of length BlockLength, one matrix per block

        % Start with the first observation in the overlapping period.
        while i<DifferentInd(end) % as long as you are in the time interval in which two individuals are tagged

            % 1. Check that the observation after the one being observed, does not fall off the end of the data set.
            i1 = 1;                                             % Index for the matrix of Start indices
            while (i+1)<size(DOld,1)
             
                % 2. Only keep blocks whose first observation is within the first half of the block
                if BlockStart > 12 % If the BlockStart is after midday
                    BlockHalf = BlockStart + BlockLength*24/2 - 24; % then the second half of the block is after midnight of the next day
                    while (DOld(i,2) > BlockHalf*60*60) && (DOld(i,2) < BlockHalf*60*60 + BlockLength*24*60*60/2) && ((i+1) < size(DOld,1)) % while the beginning of the block is in the second half of the block, but is not falling off the end of the dataset
                        i = i+1; % increase i and therefore discard previous observations
                    end
                else
                    BlockHalf = BlockStart + BlockLength*24/2; % otherwise the second half of the block is before midnight on the same day
                    while ((DOld(i,2) > BlockHalf*60*60) && ((i+1) < size(DOld,1))) || ((DOld(i,2) < BlockHalf*60*60 + BlockLength*24*60*60/2 - 24*60*60) && ((i+1) < size(DOld,1))) % while the beginning of the block is in the second half of the block
                        i = i+1; % increase i and therefore discard previous observations
                    end
                end
                if (i+1) < size(DOld,1)
                    BlockIndicesTemp(i1,1) = i;
                    BlockStartDay = DOld(i,1);
                    i1 = i1+1;
                    i = i+1;
                    while DOld(i,1) < BlockStartDay+BlockLength && i<size(DOld,1)% Move to the first observation of the day within which the next block will start
                        i = i+1; % If BlockStart does not start at the fist observation of the next day, then the next BlockStart is found on line 62 - 64
                    end
                else
                    break
                end
            end
                
            % 3. Once the start indeces of the blocks have been found, find the last observations of the blocks.
            BlockIndicesTemp(1:(end-1),2) = BlockIndicesTemp(2:end,1)-1;
            BlockIndicesTemp(end,2) = size(DOld,1);
            i3 = 1;
            for i2 = 1:size(BlockIndicesTemp,1)

                % 4. Only keep the block if the last observation of the block is after half way through the block
                if BlockStart > 12 % If the BlockStart is after midday
                    BlockHalf = BlockStart + BlockLength*24/2 - 24;
                    if (DOld(BlockIndicesTemp(i2,2),2)>BlockHalf*60*60) && (DOld(BlockIndicesTemp(i2,2),2)<BlockStart*60*60) % if the end of the block is in the second half of the block
                        BlockStartLength(i3,1) = BlockIndicesTemp(i2,1);
                        BlockStartLength(i3,2) = BlockIndicesTemp(i2,2)-BlockIndicesTemp(i2,1)+1;
                        i3 = i3+1;
                    end
                else
                    BlockHalf = BlockStart + BlockLength*24/2;
                    if ((DOld(BlockIndicesTemp(i2,2),2) > BlockHalf*60*60) || (DOld(BlockIndicesTemp(i2,2),2) < BlockStart*60*60)) % if the beginning of the block is in the second half of the block
                        BlockStartLength(i3,1) = BlockIndicesTemp(i2,1);
                        BlockStartLength(i3,2) = BlockIndicesTemp(i2,2)-BlockIndicesTemp(i2,1)+1;
                        i3 = i3+1;
                    end
                end
            end
            
            % 5. Only keep the blocks for which both individuals are tagged throughout the block within a rolling period of length Threshold
            Save = 1; % This will change to 0 if one of the intervals has less than both individuals tagged
            i5   = 1;
            for i4 = 1:size(BlockStartLength,1)
                j = BlockStartLength(i4,1);
                while j < BlockStartLength(i4,1)+BlockStartLength(i4,2)
                    l = 1;
                    while ((DOld(j+l,6)-DOld(j,6))<Threshold*60*60) && (j+l<size(DOld,1)) % Check how many observations are in the interval of at least length Threshold and make sure we do not fall off the end of the dataset
                        l = l+1;
                    end
                    if numel(unique(DOld(j:(j+l-1),5))) < 2 % if less than both individuals are tagged, -1 because the last observation is more than Threshold away from j
                        Save = 0;
                        break
                    end
                    j = j+1;
                end
                if Save == 1
                    BlockIndices(i5,1) = BlockStartLength(i4,1);
                    BlockIndices(i5,2) = BlockStartLength(i4,1)+BlockStartLength(i4,2)-1;
                    i5 = i5+1;
                end
                Save = 1;
            end
            
            % 6. otherwise save in D (matrix of all observations being used) and Divided (array with a matrix per block)
            Divided = NaN(1,6,size(BlockIndices,1));
            for i5 = 1:size(BlockIndices,1)
                D((size(D,1)+1):(size(D,1)+BlockIndices(i5,2)-BlockIndices(i5,1))+1, :) = DOld(BlockIndices(i5,1):BlockIndices(i5,2), :);
                Divided(1:(BlockIndices(i5,2)-BlockIndices(i5,1)+1), :, i5) = DOld(BlockIndices(i5,1):BlockIndices(i5,2), :);
            end
        end
    end
    