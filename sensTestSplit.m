%% Sensitivity data separation
function [sensTime,sensData,sensConcentration,sensIteration,...
    sensInterval,sensLength,sensAnalyte,sensHour] = sensTestSplit(TimePoint,...
    Measurement,testConcentration,testIteration,testInterval,...
    testLength,testAnalyte,testType,testHour)

%This function removes non-sensitivity data from the import data and
%returns new variables

sensTime = zeros(sum(strcmp(testType,'Sensitivity')),size(TimePoint,2));
sensData = zeros(sum(strcmp(testType,'Sensitivity')),size(TimePoint,2));

sensConcentration = zeros(sum(strcmp(testType,'Sensitivity')),1);
sensIteration = zeros(sum(strcmp(testType,'Sensitivity')),1);
sensInterval = zeros(sum(strcmp(testType,'Sensitivity')),1);
sensLength = zeros(sum(strcmp(testType,'Sensitivity')),1);
sensAnalyte = cell(sum(strcmp(testType,'Sensitivity')),1);
sensHour = zeros(sum(strcmp(testHour,'Sensitivity')),1);

sensCount = 1;

for count=1:1:size(Measurement,1)
    
    % Process sensitivity
    if strcmp(testType{count},'Sensitivity') == 1
        
        sensTime(sensCount,:) = TimePoint(count,:);
        sensData(sensCount,:) = Measurement(count,:);
        
        sensConcentration(sensCount,1) = testConcentration(count);
        sensIteration(sensCount,1) = testIteration(count);
        sensInterval(sensCount,1) = testInterval(count);
        sensLength(sensCount,1) = testLength(count);
        sensAnalyte{sensCount,1} = testAnalyte{count};
        sensHour(sensCount,1) = testHour(count);
        
        sensCount = sensCount + 1;
    end

end

end