%% This function checks where there are gaps of more than 12 hours and if
%  there is a gap of more than 12 hours, how long the gap is.

function [DeltaTime, GapLengths, Gaps] = GapFun(D, Gap)
    
    Time1     = D(1:(end-1), size(D,2)); % the continuous date/time counter is in the last column
    Time2     = D(2:end, size(D,2));
    DeltaTime = (Time2-Time1)/(60*60);

    Gaps(3,1)                  = D(1,1); % The first column of Gaps will just indicate the start date and time
    Gaps(4,1)                  = D(1,2)/(60*60);
    Idx                        = 1:length(DeltaTime);
    Temp                       = Idx(DeltaTime>=Gap);
    Gaps(1,2:(size(Temp,2)+1)) = Temp; % Save the index of the gaps - first row of Gaps
    Gaps(2,2:(size(Temp,2)+1)) = DeltaTime(DeltaTime>=Gap); % Saves the length of the gap in hours - second row of Gaps
    DTemp                      = D(1:(end-1), 1);
    Gaps(3,2:(size(Temp,2)+1)) = DTemp(DeltaTime>=Gap)'; % Saves the day the gap starts on - third row of Gaps
    DTemp                      = D(1:(end-1), 2)/(60*60);
    Gaps(4,2:(size(Temp,2)+1)) = DTemp(DeltaTime>=Gap)'; % and the time in hours just before the gap starts - fourth row of Gaps
    GapsTemp                   = Gaps(3:4,2:end) - Gaps(3:4,1:(end-1));
    Gaps(5:6, 2:end)           = GapsTemp; % Save the number of days and hours passed since the last gap - fifth and sixth row of Gaps
    
    % How often do which gaps appear
    GapLengths = unique(round(Gaps(2,2:(size(Temp,2)+1))));
    for i=1:length(GapLengths)
        GapLengths(2,i) = sum(round(Gaps(2,:))==GapLengths(1,i));   % counts how often gap i appears
    end
    
    GapLengths(3,1) = (D(end,6)-D(1,6))/(60*60*Gap); % Calculate how many Gaps (12 hours) the data set is long
    GapLengths(3,2) = sum(GapLengths(1,:).*GapLengths(2,:))/Gap; % how many Gaps (12 hours) are the Gaps long in total
    GapLengths(3,3) = GapLengths(3,2)/GapLengths(3,1); % Proportion of Gaps to total data set
        
end
