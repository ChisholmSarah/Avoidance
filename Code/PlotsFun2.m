%% This function plots the distances vs the significance per dyad, the
% location of the significant distances with an estimation of their home
% range, and the distance vs time for each dyad.

function PlotsFun2(Name1, Name2, DistpLess, DistpMore, DistLinInt, DataArrayLinInt, SigLevel)

    DataMatrix1 = DataArrayLinInt(DataArrayLinInt(:,5,1)~=0,:,1);
    DataMatrix2 = DataArrayLinInt(DataArrayLinInt(:,5,2)~=0,:,2);
    
    if ~isempty(DataMatrix1)
        %% P-value plot

        % Create figure name
        FigureName1  = sprintf('P-Values for %s and %s', Name1, Name2);

        % Duplicate p-values so that the lines in the plots are horizontal
        xValue1(1,1) = 0;
        xValue1(2,1) = DistpLess(1,1);
        yValue1(1,1) = DistpLess(1,2);
        yValue1(2,1) = DistpLess(1,2);
        xValue2(1,1) = 0;
        xValue2(2,1) = DistpMore(1,1);
        yValue2(1,1) = DistpMore(1,2);
        yValue2(2,1) = DistpMore(1,2);
        for ii=2:(length(DistpLess(:,1))-1) % Because Inf is the last Dist
            xValue1(2*ii-1,1) = DistpLess(ii-1,1);
            xValue1(2*ii,1)   = DistpLess(ii,1);
            yValue1(2*ii-1,1) = DistpLess(ii,2);
            yValue1(2*ii,1)   = DistpLess(ii,2);
            xValue2(2*ii-1,1) = DistpMore(ii-1,1);
            xValue2(2*ii,1)   = DistpMore(ii,1);
            yValue2(2*ii-1,1) = DistpMore(ii,2);
            yValue2(2*ii,1)   = DistpMore(ii,2);
        end

        % Create the figure
        figure('Name', FigureName1);
        hold on
        plot(xValue1, yValue1, 'r*-');
        plot(xValue2, yValue2, 'bo-');
        plot(xValue1, ones(length(xValue1))*(SigLevel/size(DistpLess,1)), 'k--');
        xlabel('Distance in meters');
        ylabel('p-value');
        title(FigureName1);
        legend('Less often', 'More often');
        hold off

        %% Time series plot
        FigureName2  = sprintf('Time Series Plot for %s and %s', Name1, Name2);

        % Create the figure
        figure('Name', FigureName2);
        plot(DataMatrix1(:,1), DistLinInt);
        datetick('x', 'dd/mm/yy');
        xlim([min(DataMatrix1(:,1)), max(DataMatrix1(:,1))]);
        xlabel('Date');
        ylabel('Distance in Meters');
        title(FigureName2);

        %% For Locations that are more often in the distances then expected
        if length(DistLinInt)>2 % don't draw location matrices if there are only two location measurements, because then all p-values are 0 anyway
            prop = 0.02; % shrink the map as well as the locations plotted
            Rows = 1:size(DistpLess,1);
            SigDistMore = Rows(~isnan(DistpMore(:,3)));
            SigDistLess = Rows(~isnan(DistpLess(:,3)));

            for i=SigDistMore
                % Create figure name
                if i==1
                    SigDist1    = 0;
                else
                    SigDist1    = DistpMore(i-1,1);
                end
                SigDist2     = DistpMore(i,1);
                FigureName3  = sprintf('More Often: Locations when %s and %s are within %d to %d meters of each other', Name1, Name2, SigDist1, SigDist2);

                % Create the figure
                xLocations      = [];
                yLocations      = [];
                xLocations(:,1) = DataMatrix1(DistLinInt>SigDist1 & DistLinInt<SigDist2,3);
                xLocations(:,2) = DataMatrix2(DistLinInt>SigDist1 & DistLinInt<SigDist2,3);
                yLocations(:,1) = DataMatrix1(DistLinInt>SigDist1 & DistLinInt<SigDist2,4);
                yLocations(:,2) = DataMatrix2(DistLinInt>SigDist1 & DistLinInt<SigDist2,4);
                ConvHull1       = convhull(DataMatrix1(:,3), DataMatrix1(:,4));
                ConvHull2       = convhull(DataMatrix2(:,3), DataMatrix2(:,4));

                figure('Name', FigureName3);
                hold on
                scatter(xLocations(:,1)*prop, yLocations(:,1)*prop, 'r*');
                scatter(xLocations(:,2)*prop, yLocations(:,2)*prop, 'bo');
                plot(DataMatrix1(ConvHull1,3)*prop, DataMatrix1(ConvHull1,4)*prop, 'r');
                plot(DataMatrix2(ConvHull2,3)*prop, DataMatrix2(ConvHull2,4)*prop, 'b');
                xlabel('meters')
                ylabel('meters')
                title(FigureName3);
                hold off

            end
        end

        %% For Locations that are more often in the distances then expected

        for i = SigDistLess
            % Create figure name
            if i == 1
                SigDist1    = 0;
            else
                SigDist1    = DistpLess(i-1,1);
            end
            SigDist2     = DistpLess(i,1);
            FigureName4  = sprintf('Less Often: Location when %s and %s are within %d to %d meters of each other', Name1, Name2, SigDist1, SigDist2);

            % Create the figure
            xLocations      = [];
            yLocations      = [];
            xLocations(:,1) = DataMatrix1(DistLinInt>SigDist1 & DistLinInt<SigDist2,3);
            xLocations(:,2) = DataMatrix2(DistLinInt>SigDist1 & DistLinInt<SigDist2,3);
            yLocations(:,1) = DataMatrix1(DistLinInt>SigDist1 & DistLinInt<SigDist2,4);
            yLocations(:,2) = DataMatrix2(DistLinInt>SigDist1 & DistLinInt<SigDist2,4);
            ConvHull1       = convhull(DataMatrix1(:,3), DataMatrix1(:,4));
            ConvHull2       = convhull(DataMatrix2(:,3), DataMatrix2(:,4));

            figure('Name', FigureName4);
            hold on
            scatter(xLocations(:,1)*prop, yLocations(:,1)*prop, 'r*');
            scatter(xLocations(:,2)*prop, yLocations(:,2)*prop, 'bo');
            plot(DataMatrix1(ConvHull1,3)*prop, DataMatrix1(ConvHull1,4)*prop, 'r');
            plot(DataMatrix2(ConvHull2,3)*prop, DataMatrix2(ConvHull2,4)*prop, 'b');
            title(FigureName4);
            hold off

        end
    end
    
end
