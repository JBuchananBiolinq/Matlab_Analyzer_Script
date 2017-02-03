%% Selectivity Plotting Function
function selPlotting(graphIndex,sensorID,TimePoint,Measurement,...
    concentration,testLength,testInterval,testIteration,testAnalyte,testInterferent,testHour)

global V summedInterferenceAverage summedInterferenceLabel summedInterferenceHour

% Process measurement data if present
if size(Measurement,1) > 1
    
    hourIndex = unique(testHour);
    % If multiple calibration/test times, split data by hour/iteration
    if length(hourIndex) > 1 && isnan(hourIndex(1)) == 0
        % Modification to compensate for different sampling rates - 11/16/15
        for i=1:1:size(Measurement,1)
            if V == 1
                splitHourInterferent(i) = Measurement(i,5/testInterval(i));
            else
                splitHourInterferent(i) = Measurement(i,size(Measurement,2));
            end
        end
        
        splitHourInterferent = splitVectorCat(splitHourInterferent,testHour);
        splitHourIndex = splitCellNum(testInterferent,testHour);
        splitHour = zeros(size(splitHourInterferent));
        for i=1:1:length(hourIndex)
            splitHour(:,i) = hourIndex(i);
        end
    % If only one calibration/test times, use raw measurement data
    else
        splitHourInterferent = Measurement(:,length(Measurement));
        splitHourIndex = testInterferent;
        splitHour = zeros(size(testInterferent));
        for i=1:1:length(hourIndex)
            splitHour(:,i) = hourIndex(i);
        end
    end
    
    labels = cell(length(unique(testInterferent)),1);
    
    % Cycle through all calibration times/instances, extract average and
    % deviation for interference runs
    loopStart = 1;
    for j=1:1:size(splitHourInterferent,2)
        
        % If no test hour is present, extract data
        if isnan(unique(testHour)) == 0
            splitHourSingle = splitHourInterferent(:,j);
            splitHourIndexSingle = splitHourIndex(:,j);
            splitHourSingle(isnan(splitHourSingle) == 1) = [];
            splitHourIndexSingle(cellfun(@isempty,splitHourIndexSingle)) = [];
            splitInterferent = splitVectorCat(splitHourSingle,splitHourIndexSingle);
        % If test hour value is present, split previous vectors by
        % interferent label
        else
            if size(splitHourInterferent,2) > 1
                splitInterferent = splitVectorCat(splitHourInterferent(:,length(splitHourInterferent)),testInterferent);
            else
                splitInterferent = splitVectorCat(splitHourInterferent(:,size(splitHourInterferent,2)),testInterferent);
                splitHourIndexSingle = splitHourIndex;
            end
            
            splitHourIndexSingle(cellfun(@isempty,splitHourIndex)) = [];
        end
        
        for i=1:1:size(splitInterferent,2)
            average(loopStart) = nanmean(splitInterferent(:,i));
            deviation(loopStart) = nanstd(splitInterferent(:,i));
            hour(loopStart) = hourIndex(j);
            loopStart = loopStart + 1;
        end
        
        if j == 1
            labels = unique(testInterferent);
        else
            labels = cat(1,labels,unique(testInterferent));
        end
        
        if size(splitHourInterferent,2) > 1
            index(j) = find(strcmp(unique(splitHourIndexSingle),'Control') == 1);
        else
            index = find(strcmp(unique(splitHourIndexSingle),'Control') == 1);
        end
        
        if isempty(index)
            close all;
            error('No Selectivity files labeled as "Control" (case sensitive)');
        end
        
    end
      
    % Plot interference data on bar plot with red error bars
    fig = figure(graphIndex);
    subplot(2,4,7); hold on; grid on;
    set(fig,'Position', [200, 200, 1440, 740]);
    colorList1 = {'k','b','g','c','m','y','r'}; %bar
    
    labelCount = length(unique(labels));
    baseline = abs(average(index(1)));
    for i=1:1:size(splitHourInterferent,2)
        bar((i-1)*labelCount+1:i*labelCount,100*average((i-1)*labelCount+1:i*labelCount)/baseline,colorList1{i});
    end
    
    h = errorbar(100*average/baseline,100*deviation/baseline,'.');
    set(h(1),'color','r');
    xlim([0 length(labels)+1]);
    set(h(1),'LineStyle','none','LineWidth',2);
    title(['Selectivity Plot for Sensor ',num2str(sensorID)]);
    ylabel('Percent Current of Initial Control');
    labels = {'',labels{:},''};
    set(gca,'xticklabel',labels)
    set(gca, 'FontSize', 10)
    %xticklabel_rotate(1:length(labels),90,labels);
    
    % Add electrode interference data to global vector
    summedInterferenceHour = vertcat(summedInterferenceHour,hour.');
    summedInterferenceAverage = vertcat(summedInterferenceAverage,average.');
    summedInterferenceLabel = vertcat(summedInterferenceLabel,labels);
    
% If no selectivity data is available, display blank plot
else 
    figure(graphIndex);
    ax = subplot(2,4,7);
    text(0.5,0.5,'No Selectivity Data',...
        'HorizontalAlignment', 'center', 'FontSize', 14);
    set(ax,'visible','off');
end

end