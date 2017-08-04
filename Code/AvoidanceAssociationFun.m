%% This function tests whether two individuals avoid each
% other or associate more often with each other than expected by
% chance.

% INPUT:    Id1     - Name of one of the individuals in single quotation
%                     marks (for file naming purposes). For example 'John'.
%           Id2     - Name of the other individual (for file naming
%                     purposes). 
%           Data1   - a 4xn1 cell array, where n1 is the number of
%                     observations for individual 1. This matrix should
%                     contain the location data for the first individual.
%                     The first column should be the date ('dd/mm/yyyy' -
%                     the single quotation marks are important), the second
%                     column should be the time ('HH:MM:SS' - the single
%                     quotation marks are important), the third and fourth
%                     column are the x and y location in meters from a
%                     reference point (which should be the same point for
%                     both individuals).
%           Data2   - a 4xn2 matrix, where n2 is the number of observations
%                     for individual 2. The matrix should have the same
%                     format as Data1.
%           Dist    - intervals to be examined (Do not include 0 and Inf,
%                     the interavals will be 0-Dist(1,1);
%                     Dist(1,1)-Dist(1,2); ...; Dist(1,end)-Inf; For
%                     example input [50, 100, 150, 200] will test 0-50m,
%                     50m-100m, 100m-150m, 150m-200m and 200m-Inf.
%           varargin: Use name value pairs, for example if you want to
%                   set SigLevel to 0.01, include "'SigLevel', 0.01"
%                   after including Id1, Id2, Data1, Data2 and Dist.
%               SigLevel    - Significance value to be used for testing.
%                             0.05 is the default value. The test adjusts
%                             for multiple testing, so this does not need
%                             to be taken into acount when choosing the
%                             significance level.
%               BlockStart  - Time of day the block will start in hours
%                             since midnight (for nocturnal animals this
%                             could be midday, as individuals are least
%                             active during that time). By default this
%                             will be midday.
%               perm        - number of permutations to be performed. 
%                             Default is 10,000.
% OUTPUT:   DistpLess   - A matrix with 3 columns. The first column
%                         gives the distances (Dist) tested. The
%                         second column gives the p-values related to
%                         the hypothesis that the two individuals are
%                         less often with the tested distances of each
%                         other. The third column gives only the
%                         p-values for the significant results and NaN
%                         for the not significant results.
%           DistpMore   - Similarly to DistpLess, but this matrix
%                         gives the p-values related to the hypothesis
%                         that the individuals are more often within
%                         the tested distances of each other.
%           p-vlaue plots              - Illustration of the p-values
%                                        for the distances (Dist)
%                                        tested, related to the
%                                        hypothesis that the
%                                        individuals were less often
%                                        and more often within the
%                                        distances of each other than
%                                        expected by chance.
%           Time series plots          - Time series of the distance
%                                        between the 2 individuals.
%           Significant location plots - Locations when the two
%                                        individuals were within a
%                                        significant distance of each
%                                        other. The convex hull of the
%                                        locations of the two
%                                        individuals show a rough
%                                        estimate of the individuals'
%                                        territories.

function [DistpLess, DistpMore] = AvoidanceAssociationFun(Id1, Id2, Data1, Data2, Dist, varargin)

    %% Set default values and replace if the values are given in varargin
    SigLevel    = 0.05;
    BlockLength = 1; % in days
    BlockStart  = 12; % in hours
    perm        = 10000;
    
    for k = 1:2:length(varargin)
        if strcmpi(varargin{k}, 'SigLevel')
            SigLevel = varargin{k+1};
        elseif strcmpi(varargin{k}, 'BlockStart')
            BlockStart = varargin{k+1};
        elseif strcmpi(varargin{k}, 'perm')
            perm = varargin{k+1};
        end
   end
        
    %% Change dates into datenum and times into seconds since midnight and save all in the variable Data
    Data = zeros(max(size(Data1,1), size(Data2,1)), size(Data1,2)-4, 2); % Number of columns is size(Data1,2)-4, because the original matrix has 3 columns for the date (YYYY, MM, DD) and 3 for the time (HH, MM, SS), but the Temp matrix only has one column for the dates and one for the times
    Data(1:size(Data1,1),1,1) = datenum(Data1(:,1:3));
    Data(1:size(Data2,1),1,2) = datenum(Data2(:,1:3));
    Data(1:size(Data1,1),2,1) = Data1(:,4)*60*60+Data1(:,5)*60+Data1(:,6); % Convert columns of HH, MM, SS to number of seconds since midnight
    Data(1:size(Data2,1),2,2) = Data2(:,4)*60*60+Data2(:,5)*60+Data2(:,6);
    
    Data(1:size(Data1,1),3,1) = Data1(:,7);
    Data(1:size(Data2,1),3,2) = Data2(:,7);
    Data(1:size(Data1,1),4,1) = Data1(:,8);
    Data(1:size(Data2,1),4,2) = Data2(:,8);
    Data(1:size(Data1,1),5,1) = ones(size(Data1,1),1);
    Data(1:size(Data2,1),5,2) = 2*ones(size(Data2,1),1);

    % Sort each of the matrices by the date and time column
    Data(1:size(Data1,1),:,1) = sortrows(Data(1:size(Data1,1),:,1),[1,2]);
    Data(1:size(Data2,1),:,2) = sortrows(Data(1:size(Data2,1),:,2),[1,2]);
    FirstDay                  = min(Data(1,1,:));
    Data(1:size(Data1,1),6,1) = (Data(1:size(Data1,1),1,1)-FirstDay)*24*60*60+Data(1:size(Data1,1),2,1); % Create continuous date/time counter in 6th column
    Data(1:size(Data2,1),6,2) = (Data(1:size(Data2,1),1,2)-FirstDay)*24*60*60+Data(1:size(Data2,1),2,2);
    
    clearvars -except BlockLength BlockStart Data Dist Id1 Id2 perm SigLevel varargin

    %% Checks where the gaps are greater or equal to 12 (GapLength)
    % hours. The first row of GapLengths are the length of the gaps
    % that appeared in the data, the second row is how often those
    % gaps appeared and the third row gives how often GapLengths fit
    % into the data set, how many GapLengths fit into the Gaps in
    % total and the proportion of the total length of the Gaps to the
    % total length of the data set.

    GapLength = 12; % A gap is if two position fixes are over 12 hours apart.
    
    % For each individual
    for jGap = 1:2
        DTemp                  = Data(:,:,jGap);
        DTemp2                 = DTemp(DTemp(:,size(Data,2)-1)~=0, :); % Only keep the non 0 entries
        [~, GapLengthsTemp, ~] = GapFun(DTemp2, GapLength);
        GapLengths(1:size(GapLengthsTemp,1), 1:size(GapLengthsTemp,2), jGap) = GapLengthsTemp;
    end

    clearvars -except BlockLength BlockStart Data Dist Id1 Id2 perm SigLevel varargin

    %% Combines the data of two individuals and divides it into
    % intervals (BlockLength). If the first measurements of an
    % interval are before BlockStart they are discarded. The interval
    % runs from BlockStart to BlockStart and for the permutations it
    % is important that only full intervals are used. The same is done
    % if the last measurements are after BlockStart.

    Threshold = 13;

    % For each pair of individuals
    DTemp = Data(:,:,1);
    D1    = DTemp(DTemp(:, size(Data,2)-1)~=0, :); % Keep only non 0 entries
    DTemp = Data(:,:,2);
    D2    = DTemp(DTemp(:, size(Data,2)-1)~=0, :);

    [~, ~, ~, ~, Divided, ~, ~] = PairFun(D1, D2, Id1, Id2, Threshold, BlockStart, BlockLength);

    clearvars -except BlockStart Dist Divided Id1 Id2 perm SigLevel varargin

    %% Linear Interpolation is used to find the expected location of each
    % individual at each time point.

    [LinIntData, ~, DataArrayLinInt] = LinIntFun(Divided, 15);

    %% Calculate the distance between the two individuals

    if ~isempty(DataArrayLinInt)
        % Calculate the distance between the two individuals
        [DistLinInt, DistCountLinInt] = DistanceFun(DataArrayLinInt, Dist);
        % Permute the blocks and calculate the distances between
        % the permuted individuals
%        [~, pLess, pMore, ~, ~, ~, ~, ~, ~] = PermutationFun(BlockStart, LinIntData, DataArrayLinInt, perm, Dist, DistCountLinInt);
        [~, pLess, pMore, ~, ~, ~, ~, ~, ~] = PermutationFun(LinIntData, DataArrayLinInt, perm, Dist, DistCountLinInt);
    else
        warning('%s and %s do not have a Matrix to Linearly Interpolate over', Id1, Id2)
    end

    if ~isempty(pLess)
        DistpLess(1:size(Dist,2)+1, 2) = pLess;
        DistpMore(1:size(Dist,2)+1, 2) = pMore;
    end

    % Let the first column of DistpLess and DistpMore give the
    % distances, the second the p-values and the third the
    % p-values that are less than the significance level including
    % the Bonferoni Correction
    DistpLess(1:size(Dist,2), 1) = Dist';
    DistpMore(1:size(Dist,2), 1) = Dist';
    DistpLess(end, 1)            = Inf;
    DistpMore(end, 1)            = Inf;
    DistpLess(:, 3)              = DistpLess(:, 2);
    DistpMore(:, 3)              = DistpMore(:, 2);
    SigLessDist                  = DistpLess(:,2)<=(SigLevel/(size(Dist,2)*2));
    SigMoreDist                  = DistpMore(:,2)<=(SigLevel/(size(Dist,2)*2));
    DistpLess(~SigLessDist, 3)   = NaN;
    DistpMore(~SigMoreDist, 3)   = NaN;

    clearvars -except DataArrayLinInt Dist DistLinInt DistpLess DistpMore Id1 Id2 perm SigLevel varargin

    %% Create plots of distance vs p-value, location for the
    % significant distances and a time series plot of distances

    PlotsFun2(Id1, Id2, DistpLess, DistpMore, DistLinInt, DataArrayLinInt, SigLevel);
    
end
