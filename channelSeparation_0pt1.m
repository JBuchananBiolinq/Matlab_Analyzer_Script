function [chanMeasurement,chanTimePoint,chanConcentration,chanIteration,...
    chanInterval,chanLength,chanHour,chanType,chanAnalyte,chanInterferent]...
    = channelSeparation_0pt1(channelIDlist,channelID,Measurement,TimePoint,testConcentration,...
    testIteration,testInterval,testLength,testHour,testType,testAnalyte,testInterferent)

% Initialize measurement channel variables
individualCount = 0;
channelCount = length(channelID(channelID == channelIDlist));
channelLength = size(TimePoint,2);

chanMeasurement = zeros(channelCount,channelLength);
chanTimePoint = zeros(channelCount,channelLength);

chanConcentration = zeros(channelCount,1);
chanIteration = zeros(channelCount,1);
chanInterval = zeros(channelCount,1);
chanLength = zeros(channelCount,1);
chanHour = zeros(channelCount,1);

chanType = cell(channelCount,1);
chanAnalyte = cell(channelCount,1);
chanInterferent = cell(channelCount,1);

for count=1:1:size(Measurement,1)
    
    % extract first channel from aggregated data
    if channelID(count) == channelIDlist(1)
        individualCount = individualCount + 1;
        
        chanMeasurement(individualCount,:) = Measurement(count,:);
        chanTimePoint(individualCount,:) = TimePoint(count,:);
        
        chanType{individualCount} = testType{count};
        chanInterferent{individualCount} = testInterferent{count};
        chanAnalyte{individualCount} = testAnalyte{count};
        chanConcentration(individualCount) = testConcentration(count);
        chanIteration(individualCount) = testIteration(count);
        chanInterval(individualCount) = testInterval(count);
        chanLength(individualCount) = testLength(count);
        chanHour(individualCount) = testHour(count);
    end
    
end

end