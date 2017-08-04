%% This function calculates the distance between each dyad and then counts
% how often the dyads are within Dist distances of each other. This function
% ignores all the cases in which there were not enough points to calculate
% the linear interpolation.

function [Distance, DistancesCount] = DistanceFun(DataArray, Distances)

    x1  = DataArray(DataArray(:,5,1)~=0,:,1); % choose all the non-0 rows
    x2  = DataArray(DataArray(:,5,2)~=0,:,2); % for each individual
    if size(x1)==size(x2) % This is just a check, but this should always be the case, otherwise something went wrong in LinIntFun
        Distance        = sqrt((x1(:,3)-x2(:,3)).^2+(x1(:,4)-x2(:,4)).^2); % Distance between individuals
        DistancesCount  = NaN(size(Distances,2)+1,1); % +1, because the final DistanceCount is for all distances greater than Distances(end)

        % For each distance in Distances, count how often the two
        % individuals were within that distance of each other
        for k=1:(size(Distances,2)+1)
            if k==1
                DistancesCount(k) = sum(Distance<=Distances(k));
            elseif k==(length(Distances)+1)
                DistancesCount(k) = sum(Distance>Distances(k-1));
            else
                DistancesCount(k) = sum(Distance<=Distances(k) & Distance>Distances(k-1));
            end
        end
    else
        warning('Something must have gone wrong with the linear interpolation. The two individuals have different number of observations.')
    end
    
end
