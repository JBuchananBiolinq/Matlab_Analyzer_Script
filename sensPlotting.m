%% Sensitivity Plotting Function
function [sensitivitySlope,LODc,LOL,LOQc,RSQ,P] = sensPlotting(graphIndex,sensorID,totalsensors,TimePoint,Measurement,...
    concentration,testLength,testInterval,testIteration,testAnalyte,testHour,gluc_increm_num)

% Current scaling factor (1e9 = plot data in nanoamps, etc.)
scaleFactor = 1e9;

% Randomly generated color list for high hour iteration numbers

colorListRandom=[];
for i=7:testHour(end)+1
    colorListRandom = [colorListRandom;rand(1,3)];
end

% Global variables
global clarkeReal clarkePredicted
global B V
global summedSensitivity summedConcentration summedDeviation summedHour summedSampleCurrentCSensitivity

% Initialize sensitivity readings for each time length (5,10,15,30,60s)
speed5 = zeros(size(Measurement,1),1);
speed10 = zeros(size(Measurement,1),1);
speed15 = zeros(size(Measurement,1),1);
speed30 = zeros(size(Measurement,1),1);
speedDev = zeros(size(Measurement,1),1);
speedEnd = zeros(size(Measurement,1),1);

% Proceed if measurement vector is nonzero
if size(Measurement,1) > 1
    
    % Extract sensor current information
    for count=1:1:size(Measurement,1)
        data = Measurement(count,:);
        % Error check to see if data for time length is present
        try
            if V == 1
                speed5(count) = data(5/testInterval(count));
            elseif V == 2
                speed5(count) = data(5/testInterval(count));
                speed10(count) = data(10/testInterval(count));
            elseif V == 3
                speed5(count) = data(5/testInterval(count));
                speed10(count) = data(10/testInterval(count));
                speed15(count) = data(15/testInterval(count));
            elseif V == 4
                speed5(count) = data(5/testInterval(count));
                speed10(count) = data(10/testInterval(count));
                speed15(count) = data(15/testInterval(count));
                speed30(count) = data(30/testInterval(count));
            elseif V == 5
                subsection = data((testLength(count)-2)/testInterval(count):testLength(count)/testInterval(count));
                speedDev(count) = mean(diff(subsection))*(1/testInterval(count));
            end
            % 60 second (or other) ending data value
            try
                speedEnd(count) = data(testLength(count)/testInterval(count));
            catch
                speedEnd(count) = data(length(data));
            end
        catch
            error('Incorrect time length selection. Rerun program.');
        end
    end
    
    if V == 1
        splitspeed5h = splitVectorCat(speed5,testHour);
    elseif V == 2
        splitspeed10h = splitVectorCat(speed10,testHour);
    elseif V == 3
        splitspeed15h = splitVectorCat(speed15,testHour);
    elseif V == 4
        splitspeed30h = splitVectorCat(speed30,testHour);
    elseif V == 5
        splitspeedDevh = splitVectorCat(speedDev,testHour);
    end
    
    splitspeedEndh = splitVectorCat(speedEnd,testHour);
    splitConHour = splitVectorCat(concentration,testHour);
    
    splitIterationHour = splitVectorCat(testIteration,testHour);
    hourListing = unique(testHour);
    
    n = 1; %color index indicator if hour timings > 8
    % Loop through each individual calibration/hour time
    initialReg = cell(1,length(hourListing));
    P=cell(1,length(hourListing));
    for j=1:1:length(hourListing)
        j;
        splitConHour(:,j);
        splitIterationHour(:,j);
        
        splitconcen = splitVectorCat(splitConHour(:,j),splitIterationHour(:,j));
        splitspeedEnd = splitVectorCat(splitspeedEndh(:,j),splitIterationHour(:,j));

        % Based on time selection, calculate & select average & deviation values
        if V == 1
            splitspeed5 = splitVectorCat(splitspeed5h(:,j),splitIterationHour(:,j));
            
            % Workaround for offset data vectors
            for k=2:1:size(splitconcen,2)
                for l=1:1:length(splitconcen(:,k))
                    if splitconcen(l,k) > splitconcen(l,k-1)
                        splitconcen(:,k) = circshift(splitconcen(:,k),1);
                        splitspeed5(:,k) = circshift(splitspeed5(:,k),1);
                    elseif splitconcen(l,k) < splitconcen(l,k-1)
                        splitconcen(:,k) = circshift(splitconcen(:,k),-1);
                        splitspeed5(:,k) = circshift(splitspeed5(:,k),-1);
                    end
                end
            end
            
            average5 = zeros(size(splitspeedEnd(:,1),1),1);
            deviation5 = zeros(size(splitspeedEnd(:,1),1),1);
            for i=1:1:length(splitspeed5)
                average5(i) = nanmean(splitspeed5(i,:));
                deviation5(i) = nanstd(splitspeed5(i,:));
            end
            average = average5;
            deviation = deviation5;
        % 10 second initialization & calculation
        elseif V == 2
            splitspeed10 = splitVectorCat(splitspeed10h(:,j),splitIterationHour(:,j));
            average10 = zeros(size(splitspeedEnd(1,:),1),1);
            deviation10 = zeros(size(splitspeedEnd(1,:),1),1);
            for i=1:1:length(splitspeed10)
                average10(i) = nanmean(splitspeed10(i,:));
                deviation10(i) = nanstd(splitspeed10(i,:));
            end
            average = average10;
            deviation = deviation10;
        % 15 second initialization & calculation
        elseif V == 3
            splitspeed15 = splitVectorCat(splitspeed15h(:,j),splitIterationHour(:,j));
            average15 = zeros(size(splitspeedEnd(1,:),1),1);
            deviation15 = zeros(size(splitspeedEnd(1,:),1),1);
            for i=1:1:length(splitspeed15)
                average15(i) = nanmean(splitspeed15(i,:));
                deviation15(i) = nanstd(splitspeed15(i,:));
            end
            average = average15;
            deviation = deviation15;
        % 30 second initialization & calculation
        elseif V == 4
            splitspeed30 = splitVectorCat(splitspeed30h(:,j),splitIterationHour(:,j));
            average30 = zeros(size(splitspeedEnd(1,:),1),1);
            deviation30 = zeros(size(splitspeedEnd(1,:),1),1);
            for i=1:1:length(splitspeed30)
                average30(i) = nanmean(splitspeed30(i,:));
                deviation30(i) = nanstd(splitspeed30(i,:));
            end
            average = average30;
            deviation = deviation30;
        elseif V == 5
            splitspeedDev = splitVectorCat(splitspeedDevh(:,j),splitIterationHour(:,j));
            averageDev = zeros(size(splitspeedEnd(1,:),1),1);
            deviationDev = zeros(size(splitspeedEnd(1,:),1),1);
            for i=1:1:length(splitspeedDev)
                averageDev(i) = nanmean(splitspeedDev(i,:));
                deviationDev(i) = nanstd(splitspeedDev(i,:));
            end
            average = averageDev;
            deviation = deviationDev;
        % 60 (or ending) second initialization & calculation
        else
            average60 = zeros(size(splitspeedEnd(1,:),1),1);
            deviation60 = zeros(size(splitspeedEnd(1,:),1),1);
            for i=1:1:length(splitspeedEnd)
                average60(i) = nanmean(splitspeedEnd(i,:));
                deviation60(i) = nanstd(splitspeedEnd(i,:));
            end
            average = average60;
            deviation = deviation60;
        end
        
        % Sort concentrations & values if out of order, remove NaN values
        [splitConcenSorted,sortIndex] = sort(splitconcen(:,1));
        average = scaleFactor*average(sortIndex);
        deviation = scaleFactor*deviation(sortIndex);
        splitConcenSorted(isnan(splitConcenSorted)) = [];
        average(isnan(average)) = [];
        deviation(isnan(deviation)) = [];
        
        % Plot sensitivity measurements by iteration #
        % elseif statements included for backwards compatibility
        fig2 = figure(graphIndex);
        set(fig2,'Position', [100, 200, 1800, 740]); % 8.5x11 ratio
        subplot(2,4,[5 6]);
        hold on; grid on;
        
        colorList1 = {'ko','bo','go','co','mo','yo','ro'}; %errorbar
        colorList2 = {'k--','b--','g--','c--','m--','y--','r--'}; %regression lines
        colorList3 = {'k','b','g','c','m','y','r'}; %scatter plot
        colorList4 = {'ko-','bo-','go-','co-','mo-','yo-','ro-'}; %errorbar
        
        try
            if hourListing(j) < 7
                errorbar(splitConcenSorted,average,deviation,colorList1{j});
            else
                errorbar(splitConcenSorted,average,deviation,'o','Color',colorListRandom(n,:));
            end
        catch
            error('Incorrect hour or calibration number. Check text files.');
        end
        
        % Add calibration/sensitivity data to full data set vectors
        hourList(:,1:length(average)) = hourListing(j);
        if size(average,1) > 1
            summedSensitivity = vertcat(summedSensitivity,average);
        else
            summedSensitivity = vertcat(summedSensitivity,average.');
        end
               
        summedDeviation = vertcat(summedDeviation,deviation);
        summedConcentration = vertcat(summedConcentration,splitConcenSorted);
        summedHour = vertcat(summedHour,hourList.');
        
        % Calculate regression of sensitivity information (linear or quad)
        try
            if B == 1
                [WLfit,CLfit,RSQ(j),P{j}] = quadraticFit(splitConcenSorted,average);
                [~,~,~,P2] = quadraticFit(average,splitConcenSorted);
                
                if j == 1
                    initialReg{j} = P2;
                end
            else
                [WLfit,CLfit,RSQ(j),P{j}] = linearFit(splitConcenSorted,average);
                [~,~,~,P2] = linearFit(average,splitConcenSorted);
                if j == 1
                    initialReg{j} = P2;
                end
            end
        catch
            close all;
            error(['Issue reading Sensitivity text data. '...
                'Check Sensitivity file header and sizes.']);
        end
        
        % Extract coefficients for future use
        if B == 1
            sensitivitySlope(j) = P{j}(2);
            sensitivityIntercept = P{j}(3);
            clarkeSlope = P2(2);
            clarkeIntercept = P2(3);
            clarkeQuadratic = P2(1);
        else
            sensitivitySlope(j) = P{j}(1);
            sensitivityIntercept = P{j}(2);
            clarkeSlope = P2(1);
            clarkeIntercept = P2(2);
            clarkeQuadratic = 0;
        end
        
        % Plot regression line on same plot as errorbar
        % elseif statements for backwards compatibility
        if hourListing(j) < 7
            plot(WLfit,CLfit,colorList2{j},'LineWidth',2);
        else
            plot(WLfit,CLfit,'--','LineWidth',2,'Color',colorListRandom(n,:));
        end

        % Calculate and plot residual standard error values
        % elseif statements for backwards compatibility
        subplot(2,4,2);
        hold on;
        concentrationList = unique(concentration);
        RSE = 100*deviation./average;
        
        if hourListing(j) < 7
            scatter(concentrationList,RSE,colorList3{j});
        else
            scatter(concentrationList,RSE,'MarkerEdgeColor',colorListRandom(n,:));
        end
        
        % If values over 5%, mark with solid red circle
        scatter(concentrationList(RSE>5), RSE(RSE>5), 'r.');
        grid on;
        title('Residual Standard Error');
        ylabel('RSE (%)');
        xlabel('Concentration (mM)');
        
        % MARD Plot - calculate and plot by concentration
        subplot(2,4,3);
        hold on; grid on;
        if V == 1
            expectedEndAll = polyval(initialReg{j},splitspeed5*1e9);
        elseif V == 2
            expectedEndAll = polyval(initialReg{j},splitspeed10*1e9);
        elseif V == 3
            expectedEndAll = polyval(initialReg{j},splitspeed15*1e9);
        elseif V == 4
            expectedEndAll = polyval(initialReg{j},splitspeed30*1e9);
        elseif V == 5
            expectedEndAll = polyval(initialReg{j},splitspeedDev*1e9);
        else
            expectedEndAll = polyval(initialReg{j},splitspeedEnd*1e9);
        end
        
        for i=1:1:length(concentrationList)
            mardAll(j,i) = 100*abs(nanmean(expectedEndAll(i,:))-concentrationList(i))/abs(concentrationList(i));
        end
        
        if hourListing(j) < 7
            plot(concentrationList,mardAll(j,:),colorList4{j});
        else
            plot(concentrationList,mardAll(j,:),'Color',colorListRandom(n,:));
            n = n+1; % Last graph, increase random color index by 1
        end
        
        title('MARD Per Concentration');
        ylabel('MARD (%)');
        xlabel('Concentration (mM)');
        
        %Plot to Summary plot which contains all of the calibration plots for each
        %electrode
            fig123=figure(123);
            set(fig123,'Position', [100, 200, 1800, 740]); % 8.5x11 ratio
            fh_subplot=subplot(2,4,sensorID);
            %subaxis(5,5,sensorID, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
            
            hold on; 
            grid on;

            % Plot regression line on same plot as errorbar
            % elseif statements for backwards compatibility      
            try
                if hourListing(j) < 7
                    errorbar(splitConcenSorted,average,deviation,colorList1{j});
                else
                    errorbar(splitConcenSorted,average,deviation,'o','Color',colorListRandom(n,:));
                end
            catch
                error('Incorrect hour or calibration number. Check text files.');
            end
            
            if hourListing(j) < 7
                plot(WLfit,CLfit,colorList2{j},'LineWidth',2);
            else
                plot(WLfit,CLfit,'--','LineWidth',2,'Color',colorListRandom(n,:));
            end
             
            hold off;
    
    end
    
    %Find LOL, LOQ, LOD
    LODc = zeros(1,length(hourListing));
    LOQc = zeros(1,length(hourListing));
    LOL = zeros(1,length(hourListing));

    for j=1:1:length(hourListing)
        row_index = length(hourListing) * gluc_increm_num * (sensorID-1) + (j-1) * gluc_increm_num + 1;
        y1=@(x) polyval(P{j},x)-(mean(summedSampleCurrentCSensitivity(row_index,:)) + 3*std(summedSampleCurrentCSensitivity(row_index,:)))*10^9;
        y2=@(x) polyval(P{j},x)-(mean(summedSampleCurrentCSensitivity(row_index,:)) + 10*std(summedSampleCurrentCSensitivity(row_index,:)))*10^9;
        LODc(j) = fsolve(y1,0);
        LOQc(j) = fsolve(y2,0); 
        LOL(j) = LimitOfLinearity(unique(concentration),mean(summedSampleCurrentCSensitivity(row_index:row_index+gluc_increm_num-1,:),2));
        
        j;
        P{j};
        mean(summedSampleCurrentCSensitivity(row_index,:))*10^9;
        std(summedSampleCurrentCSensitivity(row_index,:))*10^9;
        (mean(summedSampleCurrentCSensitivity(row_index,:)) + 3*std(summedSampleCurrentCSensitivity(row_index,:)))*10^9;
        (mean(summedSampleCurrentCSensitivity(row_index,:)) + 10*std(summedSampleCurrentCSensitivity(row_index,:)))*10^9;
        LODc(j) = fsolve(y1,0);
        LOQc(j) = fsolve(y2,0);      
        
    end
   
    
    figure(fig2)
    subplot(2,4,4);
    hold on; grid on;
    
    for j=1:1:size(mardAll,1)
        mardAvg(j) = nanmean(mardAll(j,2:size(mardAll,2)));
        mardDev(j) = nanstd(mardAll(j,2:size(mardAll,2)));
    end
    
    errorbar(hourListing,mardAvg,mardDev,'k--','LineWidth',2);
    
    title('Average MARD Per Calibration');
    ylabel('MARD (%)');
    xlabel('Calibration Number');
    
    % Aggregate and approximate Clarke error grid values from linear regression
    aggconcen = NaN;
    aggspeed = NaN;
    for i=1:1:size(splitconcen,2)
        aggconcen = vertcat(aggconcen,splitconcen(:,i));
        aggspeed = vertcat(aggspeed,splitspeedEnd(:,i));
    end
    
    clarkeReal = vertcat(clarkeReal,18*aggconcen);
    clarkePredicted = vertcat(clarkePredicted,18.*((clarkeQuadratic*(1e9*aggspeed).^2)...
        +(clarkeSlope*1e9*aggspeed)+clarkeIntercept));
    
    % Determine limit of detection - 3 standard deviations above zero measurement value
    %if hourListing(j) == 0

    
    %LODc = 0;
    %for i=1:1:length(WLfit)
    %    if CLfit(i) > LOD(j)
    %        LODc(j) = WLfit(i);
    %        break;
    %    end
    %end
    %end
    
    % Add legend and axis data to calibration/sensititivity plot
    subplot(2,4,[5 6]);
    labels_data = {'Initial Data','Initial Fit','Second Data','Second Fit','Third Data',...
        'Third Fit','Fourth Data','Fourth Fit','Fifth Data','Fifth Fit','Sixth Data',...
        'Sixth Fit','Seventh Data','Seventh Fit'};
    h_legend = legend(labels_data,'Location','SouthEast');
    set(h_legend,'FontSize',8);
    
    title('Concentration vs. Current');
    if V == 5
        ylabel('Change in Measured Current (nA)');
    else
        ylabel('Measured Current (nA)');
    end
    xlabel([testAnalyte{1},' Concentration (mM)']);
    
    % Add text information to center panel of plot
    ax = subplot(2,4,1);
    text(0.5,1,['Calibration Plot for ', testAnalyte{1}, ' Sensor'],...
        'HorizontalAlignment', 'center', 'FontSize', 16);
    text(0.5,0.85,['Sensor ID: ', num2str(sensorID)],...
        'HorizontalAlignment', 'center', 'FontSize', 14 , 'Color', 'b');
    text(0.5,0.65,'Initial Calibration Statistics',...
        'HorizontalAlignment', 'center', 'FontSize', 14 , 'Color', 'k');
    if LODc > 0.5 %0.5 %mM
        text(0.5,0.5,['LOL: ', num2str(LOL), ' mM'],...
            'HorizontalAlignment', 'center', 'FontSize', 12);
        text(0.5,0.4,['LOD: ', num2str(LODc), ' mM'],...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
        text(0.5,0.3,['LOQ: ', num2str(LOQc), ' mM'],...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
    else
        text(0.5,0.5,['LOL: ', num2str(LOL), ' mM'],...
            'HorizontalAlignment', 'center', 'FontSize', 12);
        text(0.5,0.4,['LOD: ', num2str(LODc), ' mM'],...
            'HorizontalAlignment', 'center', 'FontSize', 12);
        text(0.5,0.3,['LOQ: ', num2str(LOQc), ' mM'],...
            'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    if RSQ < 0.95 % R2 value
        text(0.5,0.2,['R2: ', num2str(RSQ)],...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
    else
        text(0.5,0.2,['R2: ', num2str(RSQ)],...
            'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    if sensitivitySlope < 0.1 % nA/mM
        text(0.5,0.1,['Sensitivity: ', num2str(sensitivitySlope), ' nA/mM'],...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
    else
        text(0.5,0.1,['Sensitivity: ', num2str(sensitivitySlope), ' nA/mM'],...
            'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    set(ax,'visible','off');
    
    else
    figure(graphIndex);
    ax = subplot(2,4,[5 6]);
    text(0.5,0.5,'No Selectivity Data',...
        'HorizontalAlignment', 'center', 'FontSize', 14);
    set(ax,'visible','off');
end
  
    %Set text labels for individual calibration summary graph
    figure(123)
    subplot(2,4,sensorID)
    hold on
    title({strcat('WE ',num2str(sensorID));strcat('Sensitivity: ',num2str(sensitivitySlope), ' nA/mM')})
    xlim([0 splitConcenSorted(end)])
    
    if sensorID == 1
        xlabel('Glucose Concentration (mM)')
        ylabel('Current(nA)')
        fh_legend = legend(labels_data,'Location','Southeast');
        set(fh_legend,'FontSize',6);
    end
    hold off
end
