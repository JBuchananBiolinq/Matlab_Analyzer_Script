%% Stability Plotting Function
function stabPlotting(graphIndex,sensorID,testData,testHour,...
    testLength,testInterval,testIteration)

global sensitivitySlope

% Proceed with plotting & calculation if data is present
if size(testData,1) > 1
    
    fig2 = figure(graphIndex);
    set(fig2,'Position', [200, 200, 1800, 740]); % 8.5x11 ratio
    subplot(2,4,8);
    title(['Stability Plot for Sensor ',num2str(sensorID)]);
    xlabel('Hours Elapsed');
    ylabel('Raw Current Output');
    hold on;
    grid on;
    
    speedEnd = zeros(size(testData(1,:),1),1);
    
    for count=1:1:size(testData,1)
        speedEnd(count) = testData(count,testLength(count)/testInterval(count));
    end
    
    if size(speedEnd(1,:),1) > 2
        splitspeedEnd = splitVectorCat(speedEnd,testIteration);
    else
        splitspeedEnd = speedEnd.';
    end
    
    average = zeros(size(splitspeedEnd(1,:),1),1);
    deviation = zeros(size(splitspeedEnd(1,:),1),1);
    
    for i=1:1:length(splitspeedEnd)
        average(i) = mean(splitspeedEnd(i,:));
        deviation(i) = std(splitspeedEnd(i,:));
    end
    
    splithour = splitVectorCat(testHour,testIteration);
    [splitHourSorted,sortIndex] = sort(splithour(:,1));
    errorbar(splitHourSorted,average(sortIndex),deviation(sortIndex));
    
    % Calculate drift rates and warmup period (time to maximum)
    driftDiff = diff(average(sortIndex));
    driftMax = max(average(sortIndex));
    timeMax = find(average(sortIndex) == driftMax);
    driftHour = diff(splitHourSorted);
    driftAvg = nanmean(driftDiff(timeMax:length(driftDiff)));
    timeMax = sum(driftHour(1:timeMax-1))*60;
    driftHour = nanmean(driftHour);
    driftPerHour = driftAvg*(1/driftHour);
    
    ax = subplot(2,4,1);
    text(0.5,0.2,['Warmup: ', num2str(max(average(sortIndex))), ' A at ' ,...
        num2str(timeMax), ' minutes'], 'HorizontalAlignment', 'center', 'FontSize', 12);
    text(0.5,0.1,['Avg Drift: ', num2str(driftPerHour), ' A per hour'],...
        'HorizontalAlignment', 'center', 'FontSize', 12);
    if sensitivitySlope ~= 0
        text(0.5,0.0,['Avg Drift: ', num2str(driftAvg/(sensitivitySlope*1e-9)),...
            ' mM per ' , num2str(driftHour), ' hour'],...
            'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    set(ax,'visible','off');
else
    figure(graphIndex);
    ax = subplot(2,4,8);
    text(0.5,0.5,'No Stability Data',...
        'HorizontalAlignment', 'center', 'FontSize', 14);
    set(ax,'visible','off');
end
end