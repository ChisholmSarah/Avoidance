%% This function creates a time series such that both individuals have an
% expected location at each point which is used to calculate their expected
% distance from each other. Linear interpolation is used.

function [LinIntData, SameTime, DataArrayLinInt] = LinIntFun(Divided, m)

    n               = size(Divided);
    DataArrayLinInt = [];
    s               = 1;
    
    for i=1:n(3) % for each block
        % Include an additional dimension to be able to have seperate
        % matrices per individual
        Ids   = unique(Divided(Divided(:,5,i)~=0, 5, i));
        Data1 = [];
        Data2 = [];
        for j=1:length(Ids)
            % Create a matrix per block per individual
            x = Divided(Divided(:,5,i)==Ids(j), :, i);
            % Count how many locations there are per individual per period
            SameTime(j,s) = size(x(x(:,5)~=0,:),1); % First two rows are number of locations per individual per period and the third and fourth row are the number of locations that were recorded within m minutes of each other
            Data1(1:size(x,1),1:size(x,2), j) = x;
            % Create a second matrix with the same dates and times as the
            % original
            Data2(1:size(x,1),1:size(x,2), j) = x;
            % but with the other individual's id
            Data2(1:size(x,1),5, j) = Ids(length(Ids)+1-j);
            % and with no x or y location (these will be filled in the
            % linear interpolation)
            Data2(1:size(x,1),3:4, j) = 0;
        end
        
        % Now we have two matrices for each individual, one being the
        % original (Data1(:,:,1), Data1(:,:,2)) and the other not having
        % any x or y locations (Data2(:,:,1), Data2(:,:,2), respectively,
        % because the Id was switched). Now we combine Data1(:,:,1) with
        % Data2(:,:,2) and Data1(:,:,2) with Data2(:,:,1) to calculate the
        % linear interpolation.
        
        for j=1:length(Ids)
            x               = [Data1(:,:,j); Data2(:,:,length(Ids)+1-j)];
            x               = x(x(:,5)~=0,:); % delete all the zero rows
            LinIntDataTemp  = sortrows(x,6); % sort by cumulative time column            
            Rows            = 1:size(LinIntDataTemp(LinIntDataTemp(:,5)~=0,:),1); % Get a vector of indicies
            Zeros           = Rows(LinIntDataTemp(LinIntDataTemp(:,5)~=0,3)==0); % for the zero x,y elements, but not the zero rows
            NotZeros        = Rows(LinIntDataTemp(:,3)~=0); % and the non zero ones
            if length(NotZeros)<2 || ~any(LinIntDataTemp(NotZeros,6)-LinIntDataTemp(NotZeros(1),6)) % You can't calculate linear interpolation if there are less than two locations available, or if all the available locations are at the same time point
                s = s-1;
                if s==0
                    DataArrayLinInt = [];
                else
                    DataArrayLinInt = DataArrayLinInt(1:a1,:,:);
                end
                LinIntDataTemp = [];
                break
            else
                for k=Zeros
                    if k<NotZeros(1) % if the first observations are not known
                        DeltaT1 = LinIntDataTemp(NotZeros(1),6) - LinIntDataTemp(Rows(k),6);
                        % Use the next two available observations to
                        % estimate the first one.
                        j1 = 2;
                        % Make sure that the next two available observations are not at the same time
                        while ((LinIntDataTemp(NotZeros(j1),6) - LinIntDataTemp(NotZeros(1),6))==0)
                            j1 = j1+1;
                        end
                        DeltaT2             = LinIntDataTemp(NotZeros(j1),6) - LinIntDataTemp(NotZeros(1),6); % Time difference between the following two observations
                        LinIntDataTemp(k,3) = round(LinIntDataTemp(NotZeros(1),3) - (LinIntDataTemp(NotZeros(j1),3)-LinIntDataTemp(NotZeros(1),3))/DeltaT2*DeltaT1); % Estimate of the x-location of the first not observed location (Location of the first observation - observed velocity travelled * time travelled from the non-observed location to the first observation)
                        LinIntDataTemp(k,4) = round(LinIntDataTemp(NotZeros(1),4) - (LinIntDataTemp(NotZeros(j1),4)-LinIntDataTemp(NotZeros(1),4))/DeltaT2*DeltaT1); % and the y-location
                    elseif k>NotZeros(end)
                        % if both the last entries are not known
                        j2 = size(NotZeros,2)-1;
                        while ((LinIntDataTemp(NotZeros(end),6) - LinIntDataTemp(NotZeros(j2),6))==0)
                            j2 = j2-1;
                        end
                        DeltaT1             = LinIntDataTemp(NotZeros(end),6) - LinIntDataTemp(NotZeros(j2),6);
                        DeltaT2             = LinIntDataTemp(Rows(k),6) - LinIntDataTemp(NotZeros(end),6);
                        LinIntDataTemp(k,3) = round(LinIntDataTemp(NotZeros(end),3) + (LinIntDataTemp(NotZeros(end),3) - LinIntDataTemp(NotZeros(j2),3))/DeltaT1*DeltaT2);
                        LinIntDataTemp(k,4) = round(LinIntDataTemp(NotZeros(end),4) + (LinIntDataTemp(NotZeros(end),4) - LinIntDataTemp(NotZeros(j2),4))/DeltaT1*DeltaT2);
                    else % For all other missing observations (not first and last)
                        % Get the largest index (of Rows) that is smaller than
                        % the 0's index
                        IdxLess             = max(NotZeros(NotZeros<Rows(k)));
                        % and the smallest index (of Rows) that is smaller than
                        % the 0's index
                        IdxMore             = min(NotZeros(NotZeros>Rows(k)));
                        % calculate the time differenced from the lower observation to the unknown location
                        DeltaT1             = LinIntDataTemp(Rows(k),6) - LinIntDataTemp(Rows(IdxLess),6);
                        % and from the higher observation to the lower location
                        DeltaT2             = LinIntDataTemp(Rows(IdxMore),6)- LinIntDataTemp(Rows(IdxLess),6);
                        % Calculate the x and y linear interpolation
                        LinIntDataTemp(k,3) = round(LinIntDataTemp(Rows(IdxLess),3) + (LinIntDataTemp(Rows(IdxMore),3) - LinIntDataTemp(Rows(IdxLess),3))/DeltaT2*DeltaT1);
                        LinIntDataTemp(k,4) = round(LinIntDataTemp(Rows(IdxLess),4) + (LinIntDataTemp(Rows(IdxMore),4) - LinIntDataTemp(Rows(IdxLess),4))/DeltaT2*DeltaT1);
                    end
                end
           end

            %% Delete all the double entries (entries where both individuals'
            % location was recorded within m minutes of each other)
            % and count the number of times they were recorded at
            % exactly the same time

            xTemp    = LinIntDataTemp(LinIntDataTemp(:,5)~=0,:); % Only keep non-zero rows
            x1       = xTemp(1:(end-1),6);
            x2       = xTemp(2:end,6);
            a        = abs(x1-x2)>=(m*60); % a vector of 1s and 0s to show whether the consecutive tags were outside of m minutes of each other
            a(end+1) = 1;
            LinIntData(1:size(LinIntDataTemp(a,:),1), :,j,s) = LinIntDataTemp(a,:);
            SameTime(j+length(Ids),s) = sum(~a); % number of lacations that were measured at exactly the same time

            %% Create one matrix per individual, not per block

            if s==1 % s is the number of matrices per individual per block
                a1 = 0; % a1 is the starting matrix for this block in DataMatrixLinInt
            else
                a1 = size(DataArrayLinInt(DataArrayLinInt(:,5,j)~=0,:,j),1);
            end
            x                                   = LinIntData(LinIntData(:,5,j,s)~=0,:,j,s);
            b1                                  = size(x,1);
            DataArrayLinInt((a1+1):(a1+b1),:,j) = x;
        end
        s=s+1;
    end
    
    SameTime(size(SameTime,1)+1,1) = sum(SameTime(1,:));
    SameTime(size(SameTime,1),2)   = sum(SameTime(2,:));
    SameTime(size(SameTime,1),3)   = sum(SameTime(3,:));
    
end
