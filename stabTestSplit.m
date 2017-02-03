%% Stability data separation
function [stabTime,stabData,stabConcentration,stabIteration,...
    stabInterval,stabLength,stabAnalyte,stabHour] =...
    stabTestSplit(TimePoint,Measurement,testConcentration,testIteration,...
    testInterval,testLength,testAnalyte,testType,testHour)

stabTime = zeros(sum(strcmp(testType,'Stability')),size(TimePoint,2));
stabData = zeros(sum(strcmp(testType,'Stability')),size(TimePoint,2));

stabConcentration = zeros(sum(strcmp(testType,'Stability')),1);
stabIteration = zeros(sum(strcmp(testType,'Stability')),1);
stabInterval = zeros(sum(strcmp(testType,'Stability')),1);
stabLength = zeros(sum(strcmp(testType,'Stability')),1);
stabHour = zeros(sum(strcmp(testType,'Stability')),1);
stabAnalyte = cell(sum(strcmp(testType,'Stability')),1);

stabCount = 1;

for count=1:1:size(Measurement,1)
    
    % Process Stability
    if strcmp(testType{count},'Stability') == 1
        
        stabTime(stabCount,:) = TimePoint(count,:);
        stabData(stabCount,:) = Measurement(count,:);
        
        stabConcentration(stabCount,1) = testConcentration(count);
        stabIteration(stabCount,1) = testIteration(count);
        stabInterval(stabCount,1) = testInterval(count);
        stabLength(stabCount,1) = testLength(count);
        stabAnalyte{stabCount,1} = testAnalyte{count};
        
        stabHour(stabCount,1) = testHour(count);
        
        stabCount = stabCount + 1;
        
    end
    
end

end
