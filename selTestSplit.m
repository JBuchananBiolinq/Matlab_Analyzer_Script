%% Selectivity data separation
function [selectTime,selectData,selectConcentration,selectIteration,...
    selectInterval,selectLength,selectAnalyte,selectInterferent,selectHour] = ...
    selTestSplit(TimePoint,Measurement,testConcentration,...
    testIteration,testInterval,testLength,testAnalyte,...
    testType,testInterferent,testHour)

selectTime = zeros(sum(strcmp(testType,'Selectivity')),size(TimePoint,2));
selectData = zeros(sum(strcmp(testType,'Selectivity')),size(TimePoint,2));

selectConcentration = zeros(sum(strcmp(testType,'Selectivity')),1);
selectIteration = zeros(sum(strcmp(testType,'Selectivity')),1);
selectInterval = zeros(sum(strcmp(testType,'Selectivity')),1);
selectLength = zeros(sum(strcmp(testType,'Selectivity')),1);
selectAnalyte = cell(sum(strcmp(testType,'Selectivity')),1);
selectInterferent = cell(sum(strcmp(testType,'Selectivity')),1);
selectHour = zeros(sum(strcmp(testType,'Selectivity')),1);

selectCount = 1;

for count=1:1:size(Measurement,1)
    
    % Process Selectivity
    if strcmp(testType{count},'Selectivity') == 1
        
        selectTime(selectCount,:) = TimePoint(count,:);
        selectData(selectCount,:) = Measurement(count,:);
        
        selectConcentration(selectCount,1) = testConcentration(count);
        selectIteration(selectCount,1) = testIteration(count);
        selectInterval(selectCount,1) = testInterval(count);
        selectLength(selectCount,1) = testLength(count);
        selectAnalyte{selectCount,1} = testAnalyte{count};
        
        selectInterferent{selectCount,1} = testInterferent{count};
        selectHour(selectCount,1) = testHour(count);
        
        selectCount = selectCount + 1;
        
    end
end

end
