function [TimePoint,Current,Measurement,...
    testConcentration,testIteration,testInterval,testLength,...
    channelID,testType,testAnalyte,testInterferent,testHour]...
    = dataImport_0pt1(rawCalibDir,calibrationList)

% Preallocate data vectors for text file input
TimePoint = zeros(length(calibrationList),20);
Current = zeros(length(calibrationList),20);
Measurement = zeros(length(calibrationList),20);

testConcentration = zeros(length(calibrationList),1);
testIteration = zeros(length(calibrationList),1);
testInterval = zeros(length(calibrationList),1);
testLength = zeros(length(calibrationList),1);
channelID = zeros(length(calibrationList),1);

testType = cell(length(calibrationList),1);
testAnalyte = cell(length(calibrationList),1);
testInterferent = cell(length(calibrationList),1);
testHour = zeros(length(calibrationList),1);

% Import data for all files
measCount = 1;
for count=1:1:length(calibrationList)
    
    
    fileName = calibrationList(count).name
    brokenName = regexp(fileName,'_','split');
    testFile = [rawCalibDir,'/',fileName];
    testFileID = fopen(testFile);
    dataHeader = textscan(testFileID,'%s','delimiter','\n'); % import file into cell
    dataHeader = dataHeader{1}; % break out file into cells
    fclose(testFileID);
    
    singleTimestamp = dataHeader{1};
    singleAnalyte = brokenName{2};
    singleType = brokenName{3};
    
    % Contingency for extra header information ("Ref. Electrode" in position #6
    headerCheck = regexp(dataHeader{6},':','split');
    if strcmp(headerCheck{1},'Ref. Electrode')
        singleInterval = regexp(dataHeader{11},'=','split');
        singleInterval = str2double(singleInterval{2});
        singleLength = regexp(dataHeader{12},'=','split');
        singleLength = str2double(singleLength{2});
    else
        singleInterval = regexp(dataHeader{10},'=','split');
        singleInterval = str2double(singleInterval{2});
        singleLength = regexp(dataHeader{11},'=','split');
        singleLength = str2double(singleLength{2});
    end
    
    % Fix testing type if file has lower-case lettering
    if strcmp(singleType,'sensitivity') == 1
        singleType = 'Selectivity';
    elseif strcmp(singleType,'selectivity') == 1
        singleType = 'Selectivity';
    elseif strcmp(singleType,'stability') == 1
        singleType = 'Stability';
    end
        
    singleConcentration = 1000*str2double(brokenName{4});
    
    if strcmp(singleType,'Sensitivity') == 1
        % Short sensitivity filename without hour/calibration number
        if length(brokenName) < 5
            singleIteration = str2double(strtok(brokenName{4},'.'));
            singleInterferent = NaN;
            singleHour = 0;
        elseif length(brokenName) < 6
            singleIteration = str2double(strtok(brokenName{5},'.'));
            singleInterferent = NaN;
            singleHour = 0;
        % Normal sensitivity filename
        else
            singleIteration = str2double(strtok(brokenName{6},'.'));
            singleInterferent = NaN;
            singleHour = str2double(brokenName{5});
        end
    elseif strcmp(singleType,'Selectivity') == 1
        try
            singleInterferent = brokenName{5};
            % Selectivity file with no hour/calibration value
            if length(brokenName) < 7
                singleIteration = str2double(strtok(brokenName{6},'.'));
                singleHour = NaN;
            % Normal selectivity file name with hour/calibration value
            else
                singleIteration = str2double(strtok(brokenName{7},'.'));
                singleHour = str2double(brokenName{6});
            end
        catch
            close all;
            error(['File name missing field(s). Please check: ',...
                calibrationList(count).name]);
        end
    elseif strcmp(singleType,'Stability') == 1
        % Normal stability filename
        if length(brokenName) > 4
            singleIteration = str2double(strtok(brokenName{6},'.'));
            singleInterferent = NaN;
            singleHour = str2double(brokenName{5});
        % Short stability filename (for long tests)
        else
            singleIteration = 1;
            singleIndex = str2double(strtok(brokenName{4},'.'));
            singleInterferent = NaN;
            singleTimestamp = datevec(singleTimestamp);
            if singleIndex == 1
                startHour = singleTimestamp;
            end
            singleHour = singleTimestamp-startHour;
            singleHour = (singleHour(3)*86400+singleHour(4)*3600+...
                singleHour(5)*60+singleHour(6))/3600;
        end
    else
        error(['Test type incorrectly labeled. Please check: ',...
            calibrationList(count).name]);
    end
    
    % Check to see if software has additional channel header information
    % and find first line of data (column count > 1)
    for j=9:1:length(dataHeader)
        columnHeader = regexp(dataHeader{j},',','split');
        if size(columnHeader,2) > 1
            i = j+2;
            offset = j+1;
            break;
        end
    end
    
    channelList = regexp(brokenName{1},'-','split');

    % Import data from text file for arbitrary number of channels
    while i <= length(dataHeader)
        splitLine = regexp(dataHeader{i},',','split'); 
        for j=0:1:length(splitLine)-2
            TimePoint(measCount+j,i-offset) = str2double(splitLine{1});
            Current(measCount+j,i-offset) = str2double(splitLine{2+j});
            Measurement(measCount+j,i-offset) = str2double(splitLine{2+j});
        end
        i = i + 1;
    end
    
    % Extract test information for arbitrary number of channels
  
    for i=2:1:length(columnHeader)
        %channelList
        %i-1
        %channelList{i-1}
        %measCount
        %channelID
        %channelList
        %i
        channelID(measCount) = str2double(channelList{i-1});
        testType{measCount} = singleType; 
        testAnalyte{measCount} = singleAnalyte;
        testConcentration(measCount) = singleConcentration;
        testIteration(measCount) = singleIteration;
        testInterval(measCount) = singleInterval;
        testLength(measCount) = singleLength;
        testInterferent{measCount} = singleInterferent;
        testHour(measCount) = singleHour;
        measCount = measCount+1;
    end
    
end

end