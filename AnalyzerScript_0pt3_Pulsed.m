function AnalyzerScript_0pt3_Pulsed()
% 0.3_Pulsed: Specialized import script to handle pulsed data input
% Note: Takes considerable time to import files

% Required files for operation:
% dataImport_0pt1.m
% channelSeparation_0pt1.m
% clarke.m
% splitVectorCat.m
% quadraticFit.m
% linearFit.m

clear all; close all;

% Global variable initialization & definement
global clarkeReal clarkePredicted googleDrive 
global B %V
%global summedSensitivity summedConcentration summedCounter
%global summedHour summedElectrodeNumber
%global summedInterferenceAverage summedInterferenceLabel summedInterferenceHour
summedCounter = 1;

scaleFactor = 1e9;

% Establish main file directory
googleDrive = 'C:\Users\Public\Dropbox\Analyzer Results\';

% Select sensor data folder
rawCalibDir = uigetdir([googleDrive,'Raw Text Data\']);

% Extract directory information for reference
brokenDir = regexp(rawCalibDir,'\','split');
testName = brokenDir{length(brokenDir)};

% Select pdf output folder, create if not present
outputFolder = [googleDrive,'PDF Reports\',brokenDir{length(brokenDir)-1},'\'];
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

folderList = dir(rawCalibDir);
%dirFlags = folderList.isdir;
folderList(1:2) = [];

% Choose type of regression/fit and desired sampling time
B = menu('Quadratic or Linear Fit?','Quadratic','Linear');
C = menu('Split Odd and Even Electrodes?','Yes','No');

%A = menu('Plot Clarke Grid?','Yes','No');
clarkeReal = NaN;
clarkePredicted = NaN;

%% Split channel data & individual sensor plotting - loop through each available electrode

average = zeros(1,length(folderList));
deviation = zeros(1,length(folderList));

% Loop through individual concentration folders
for k=1:1:size(folderList)
    
    % Extract list of files from directory
    subdirectory = [rawCalibDir,'\',folderList(k).name];
    concentrationList(k) = str2double(folderList(k).name);
    calibrationList = dir([subdirectory,'\*.txt']);
    
    % Import all data from files and names of files within selected folder
    [TimePoint,~,Measurement,...
        testConcentration,testIteration,testInterval,testLength,...
        channelID,testType,testAnalyte,testInterferent,testHour]...
        = dataImport_0pt1(subdirectory,calibrationList);
    
    channelIDlist = unique(channelID);
    
    % Loop through individual channels
    for i=1:1:length(channelIDlist)

        % Extract information for individual channel
        [MeasurementC,TimePointC,~,testIterationC,...
            ~,~,~,~,~,~]...
            = channelSeparation_0pt1(channelIDlist(i),channelID,Measurement,TimePoint,testConcentration,...
            testIteration,testInterval,testLength,testHour,testType,testAnalyte,testInterferent);
        %Measurment C row # = pulse #, column number = each sample per time
        %increment
        %[~,sortIndex] = sort(testIterationC(:,1));
        
        endIndex = size(MeasurementC,1);
        sortedMeasurement = zeros(length(testIteration),testLength(1)/testInterval(1));
        sortedTimePoint = zeros(length(testIteration),testLength(1)/testInterval(1));
        for j=1:1:endIndex
            sortedMeasurement(testIterationC(j),:) = MeasurementC(j,:);
            sortedTimePoint(testIterationC(j),:) = TimePointC(j,:);
        end
        average(i,k) = scaleFactor*nanmean(sortedMeasurement((endIndex-round(endIndex/10)):endIndex,20));
        deviation(i,k) = scaleFactor*nanstd(sortedMeasurement((endIndex-round(endIndex/10)):endIndex,20));
            
    end
    
    %output = figure(i*10);
    %set(output,'PaperUnits','inches','PaperPosition',[0 0 12.96 8.64],'PaperSize',[12.96 8.64]);
    %saveas(output,[outputFolder,testName,' - ',num2str(channelIDlist(i)),'.pdf'],'pdf') %Save figure
    
end

% Open and configure plotting window
fig = figure(900);
set(fig,'Position', [100, 100, 1400, 740]); % 8.5x11 ratio
subplot(2,3,4);
hold on; grid on;

% Plot individual sensor data, calculate and plot
for i=1:1:size(average,1)
    if isnan(average(i,1)) == 0
        % Select for quadratic or linear regression
        if B == 1
            [WLfit,CLfit,~,~] = quadraticFit(concentrationList,average(i,:));
        else
            [WLfit,CLfit,~,~] = linearFit(concentrationList,average(i,:));
        end
        
        if CLfit(1) > CLfit(length(CLfit)) %Label lines red if inverted
            errorbar(concentrationList,average(i,:),deviation(i,:),'ro');
            plot(WLfit,CLfit,'r');
        else
            errorbar(concentrationList,average(i,:),deviation(i,:),'ko');
            plot(WLfit,CLfit,'k');
        end
    end
end

title('Individual Sensors - Initial Calibration');
xlabel('Concentration (mM)');

% Plot group calibration/regression based on individual sensors
subplot(2,3,[5 6]);
hold on; grid on;

% Perform average and deviation calculations on full electrode set
if C == 2 % For single regression
    for i=1:1:length(concentrationList)
        pooledAverage(i) = nanmean(average(:,i));
        pooledDeviation(i) = nanmean(deviation(:,i));
    end
    
    if B == 1
        [WLfit,CLfit,~,PP] = quadraticFit(concentrationList,pooledAverage);
    else
        [WLfit,CLfit,~,PP] = linearFit(concentrationList,pooledAverage);
    end
    
    if CLfit(1) > CLfit(length(CLfit)) %Label lines red if inverted
        errorbar(concentrationList,pooledAverage,pooledDeviation,'ro');
        plot(WLfit,CLfit,'r--','LineWidth',2);
    else
        errorbar(concentrationList,pooledAverage,pooledDeviation,'ko');
        plot(WLfit,CLfit,'k--','LineWidth',2);
    end
    
    labels = {'Electrode Data','Electrode Fit'};
else % For split regressions by odd & even electrodes
    for i=2:2:size(average,1)
        averageO(i/2,:) = average(i-1,:);
        averageE(i/2,:) = average(i,:);
        deviationO(i/2,:) = deviation(i-1,:);
        deviationE(i/2,:) = deviation(i,:);
    end
    
    for i=1:1:length(concentrationList)
        pooledAverageE(i) = nanmean(averageE(:,i));
        pooledDeviationE(i) = nanmean(deviationE(:,i));
        pooledAverageO(i) = nanmean(averageO(:,i));
        pooledDeviationO(i) = nanmean(deviationO(:,i));
    end
    
    if B == 1
        [WLfitO,CLfitO,~,~] = quadraticFit(concentrationList,pooledAverageO);
        [WLfitE,CLfitE,~,~] = quadraticFit(concentrationList,pooledAverageE);
    else
        [WLfitO,CLfitO,~,~] = linearFit(concentrationList,pooledAverageO);
        [WLfitE,CLfitE,~,~] = linearFit(concentrationList,pooledAverageE);
    end
    
    if CLfitE(1) > CLfitE(length(CLfitE)) %Label lines red if inverted
        errorbar(concentrationList,pooledAverageE,pooledDeviationE,'ro');
        plot(WLfitE,CLfitE,'r--','LineWidth',2);
    else
        errorbar(concentrationList,pooledAverageE,pooledDeviationE,'bo');
        plot(WLfitE,CLfitE,'b--','LineWidth',2);
    end
    
    if CLfitO(1) > CLfitO(length(CLfitO)) %Label lines red if inverted
        errorbar(concentrationList,pooledAverageO,pooledDeviationO,'ro');
        plot(WLfitO,CLfitO,'r--','LineWidth',2);
    else
        errorbar(concentrationList,pooledAverageO,pooledDeviationO,'ko');
        plot(WLfitO,CLfitO,'k--','LineWidth',2);
    end
    
    labels = {'Even Electrodes','Even Fit','Odd Electrodes','Odd Fit'};
end

title('Group Regressions by Time');
h_legend = legend(labels,'Location','SouthEastOutside');
set(h_legend,'FontSize',8);
xlabel('Concentration (mM)');
ylabel('Sensor Current (nA)');

%{
figure(1);
hold on;
scatter(concentrationList,pooledAverage,'r');
scatter(concentrationList,pooledAverageO,'b');
scatter(concentrationList,pooledAverageE,'g');
for i=1:1:length(concentrationList)
    scatter(concentrationList,average(i,:),'k');
end
%}

% Calculate and plot residual standard error values
figure(900);
subplot(2,3,2);
hold on; grid on;
concRSE = concentrationList;%unique(concentrationList); % If concentrations duplicated, take unique vector

if C == 2 % For single regression
    RSE = 100*pooledDeviation./pooledAverage;
    scatter(concRSE,RSE,'k');
    % If values over 5%, mark with solid red dot in center
    scatter(concRSE(RSE>5), RSE(RSE>5), 'r.');
% For split electrode values
else
    RSEe = 100*pooledDeviationE./pooledAverageE;
    RSEo = 100*pooledDeviationO./pooledAverageO;
    scatter(concRSE,RSEe,'b');
    scatter(concRSE,RSEo,'k');
    % If values over 5%, mark with solid red dot in center
    scatter(concRSE(RSEe>5), RSEe(RSEe>5), 'r.');
    scatter(concRSE(RSEo>5), RSEo(RSEo>5), 'r.');
end

grid on;
title('Residual Standard Error');
ylabel('RSE (%)');
xlabel('Concentration (mM)');

% MARD Plot - calculate and plot by concentration
subplot(2,3,3);
hold on; grid on;

if C == 2
    figure(2);
    hold on; grid on;
    
    [~,~,~,PP] = quadraticFit(pooledAverage,concentrationList);
    
    for i=1:1:size(average,1)
        expectedVal(i,:) = polyval(PP,average(i,:));
        scatter(concentrationList,expectedVal(i,:),'ko');
    end
    
else
    figure(1);
    hold on; grid on;
    
    [~,~,~,PPo] = quadraticFit(pooledAverageO,concentrationList);
    [~,~,~,PPe] = quadraticFit(pooledAverageE,concentrationList);
    
    for i=1:1:size(averageO,1)
        expectedOdd(i,:) = polyval(PPo,averageO(i,:));
        scatter(concentrationList,expectedOdd(i,:),'ko');
    end
    
    for i=1:1:size(averageE,1)
        expectedEven(i,:) = polyval(PPe,averageE(i,:));
        scatter(concentrationList,expectedEven(i,:),'bo');
    end
    
    %%{
    for i=1:1:size(expectedOdd,1)
        %errorbar(concentrationList,expectedOdd(i,:),deviation(i,:),'ko');
        plot(concentrationList,expectedOdd(i,:),'ko');
    end
    
    for i=1:1:size(expectedEven,1)
        %errorbar(concentrationList,expectedEven(i,:),deviation(i,:),'bo');
        plot(concentrationList,expectedEven(i,:),'bo');
    end
    %}
    
    for i=1:1:length(concentrationList)
        for j=1:1:size(expectedOdd,1)
            mardOdd(j,i) = 100*abs(expectedOdd(j,i)-concentrationList(i))/abs(concentrationList(i));
        end
        for j=1:1:size(expectedEven,1)
            mardEven(j,i) = 100*abs(expectedEven(j,i)-concentrationList(i))/abs(concentrationList(i));
        end
    end
    
    figure(900);
    subplot(2,3,3);
    hold on; grid on;
    
    for i=1:1:size(expectedOdd,1)
        scatter(concentrationList,mardOdd(i,:),'ko');
    end
    
    for i=1:1:size(expectedEven,1)
        scatter(concentrationList,mardEven(i,:),'bo');
    end
    
end

%{
if hourListing(j) < 7
    plot(concentrationList,mardAll(j,:),colorList4{j});
else
    plot(concentrationList,mardAll(j,:),'Color',colorListRandom(n,:));
    n = n+1; % Last graph, increase random color index by 1
end
%}

title('MARD Per Concentration');
ylabel('MARD (%)');
xlabel('Concentration (mM)');

LOD =(pooledAverage(1) + pooledDeviation(1)*3);
% Determine limit of detection
%LODc = 0;
%for i=1:1:length(WLfit)
%    if CLfit(i) > LOD
%        LODc = WLfit(i);
%        break;
%    end
%end

% Add text information to center panel of plot
figure(900);
ax = subplot(2,3,1);
text(0.5,1,['Calibration Plot for ', testAnalyte{1}, ' Sensor'],...
    'HorizontalAlignment', 'center', 'FontSize', 16);
text(0.5,0.85,testName,...
    'HorizontalAlignment', 'center', 'FontSize', 14 , 'Color', 'b');
text(0.5,0.65,'Initial Calibration Statistics',...
    'HorizontalAlignment', 'center', 'FontSize', 14 , 'Color', 'k');

if LOD > 0.1 %0.5 %mM
    text(0.5,0.4,['LOD: ', num2str(LOD), ' nA'],...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
else
    text(0.5,0.4,['LOD: ', num2str(LOD), ' nA'],...
        'HorizontalAlignment', 'center', 'FontSize', 12);
end
%{
if RSQ < 0.95 % R2 value
    text(0.5,0.3,['R2: ', num2str(RSQ)],...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
else
    text(0.5,0.3,['R2: ', num2str(RSQ)],...
        'HorizontalAlignment', 'center', 'FontSize', 12);
end
%}
%{
if sensitivitySlope < 0.1 % nA/mM
    text(0.5,0.2,['Sensitivity: ', num2str(sensitivitySlope), ' nA/mM'],...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
else
    text(0.5,0.2,['Sensitivity: ', num2str(sensitivitySlope), ' nA/mM'],...
        'HorizontalAlignment', 'center', 'FontSize', 12);
end
%}
set(ax,'visible','off');

end

%}
    
%{
%% Plot data for full set of electrodes
% Color listings for plots with hours/iterations from 1 to 7
colorList1 = {'ko','bo','go','co','mo','yo','ro'}; %errorbar
colorList2 = {'k--','b--','g--','c--','m--','y--','r--'}; %regression lines

% Randomly generated color list for high hour iteration numbers
colorListRandom = [rand(1,3);rand(1,3);rand(1,3);rand(1,3);rand(1,3);rand(1,3);rand(1,3)];

% Split all sensitivity data by electrode number
splitConcElectrodeNumber = splitVectorCat(summedConcentration,summedElectrodeNumber);
splitSensElectrodeNumber = splitVectorCat(summedSensitivity,summedElectrodeNumber);
splitHourElectrodeNumber = splitVectorCat(summedHour,summedElectrodeNumber);

% Check if data needs to be split into odd and even channel groupings
if W == 1
    
    for i=1:1:size(splitConcElectrodeNumber,2)/2
        splitConcOdd(:,i) = splitConcElectrodeNumber(:,2*i-1);
        splitConcEven(:,i) = splitConcElectrodeNumber(:,2*i);
        splitSensOdd(:,i) = splitSensElectrodeNumber(:,2*i-1);
        splitSensEven(:,i) = splitSensElectrodeNumber(:,2*i);
        splitHourOdd(:,i) = splitHourElectrodeNumber(:,2*i-1);
        splitHourEven(:,i) = splitHourElectrodeNumber(:,2*i);
    end
    
    % Calculate and plot sensitivity values for odd numbered channels
    output = figure(900);
    set(output,'Position', [50, 200, 1800, 740]); % 8.5x11 ratio
    subplot(2,4,5);
    hold on; grid on;
    [averageSensitivity,deviationSensitivity,WLfit,CLfit,RSQ,sensitivitySlope,...
        calAverage,calConcentration,calDeviation,calIteration,averageConcentration] = summaryPlot(...
        splitConcOdd,splitSensOdd,splitHourOdd,colorList1,colorList2,colorListRandom);
    
    % Calculate and plot sensitivity values for even numbered channels
    % VALUES ARE NOT USED FOR DATA SUMMARY INFORMATION AS OF 6/14/16
    output = figure(902);
    set(output,'Position', [50, 200, 1800, 740]); % 8.5x11 ratio
    subplot(2,4,5);
    hold on; grid on;
    [~,~,~,~,~,~,...
        ~,~,~,~,~] = summaryPlot(...
        splitConcEven,splitSensEven,splitHourEven,colorList1,colorList2,colorListRandom);
    
% Calculate and plot sensitivity values for all channels
else 
    output = figure(900);
    set(output,'Position', [50, 200, 1800, 740]); % 8.5x11 ratio
    subplot(2,4,5);
    hold on; grid on;
    [averageSensitivity,deviationSensitivity,WLfit,CLfit,RSQ,sensitivitySlope,...
        calAverage,calConcentration,calDeviation,calIteration,averageConcentration] = summaryPlot(...
        splitConcElectrodeNumber,splitSensElectrodeNumber,...
        splitHourElectrodeNumber,colorList1,colorList2,colorListRandom);
end

LOD =(averageSensitivity(1) + deviationSensitivity(1)*3);
% Determine limit of detection
LODc = 0;
for i=1:1:length(WLfit)
    if CLfit(i) > LOD
        LODc = WLfit(i);
        break;
    end
end

% Plot interference data
figure(900);
ax = subplot(2,4,8);
hold on;
grid on;

splitInterferenceTest = splitVectorCat(summedInterferenceAverage,summedInterferenceHour);

% Extract and plot average & deviation data for each time iteration
for i=1:1:size(splitInterferenceTest,2)
    splitLabel = unique(summedInterferenceLabel);
    splitInterference = splitVectorCat(splitInterferenceTest(:,i),summedInterferenceLabel);
    
    for k=1:1:size(splitInterference,2)
        intAverage(k) = nanmean(splitInterference(:,k));
        intDeviation(k) = nanstd(splitInterference(:,k));
    end
    
    % Initialize total arrays if i = 1; otherwise concatenate vector
    if i == 1
        totalAverage = intAverage.';
        totalDeviation = intDeviation.';
        totalLabel = splitLabel;
    else
        totalAverage = vertcat(totalAverage,intAverage.');
        totalDeviation = vertcat(totalDeviation,intDeviation.');
        totalLabel = vertcat(totalLabel,splitLabel);
    end
    
end

% Find control index at time 0 and use as normalization of other data
if isempty(splitInterferenceTest) == 0
    index = find(strcmp(unique(totalLabel),'Control') == 1);
    labelCycle = length(unique(totalLabel));
    baseline = abs(totalAverage(index));
    
    if isempty(totalLabel) == 0
        bar(1:labelCycle,100*totalAverage(1:labelCycle)./baseline,'k');
    end
    
    if length(totalLabel) > labelCycle
        bar(labelCycle+1:2*labelCycle,100*totalAverage(labelCycle+1:2*labelCycle)./baseline,'b');
    end
    
    if length(totalLabel) > labelCycle*2
        bar(2*labelCycle+1:3*labelCycle,100*totalAverage(2*labelCycle+1:3*labelCycle)./baseline,'g');
    end
    
    if length(totalLabel) > labelCycle*2
        bar(3*labelCycle+1:4*labelCycle,100*totalAverage(3*labelCycle+1:4*labelCycle)./baseline,'g');
    end
    
    h = errorbar(100*totalAverage/baseline,100*totalDeviation/baseline);
    set(h(1),'color','r');
    
    % Add additional labels to interference plot
    xlim([0 length(totalLabel)+1]);
    set(h(1),'LineStyle','none','LineWidth',2);
    title('Selectivity Plot for Sensors');
    ylabel('Percent Current of Control (Initial Cal.)');
    xticklabel_rotate(1:length(totalLabel),90,totalLabel);
else
    text(0.5,0.5,'No Selectivity Data',...
        'HorizontalAlignment', 'center', 'FontSize', 14);
    set(ax,'visible','off');
end

% RSD Plot
subplot(2,4,2);
hold on;
%expected60 = (concentrationList.^2)*P(3) + concentrationList*P(2) + P(1);
compressedAvg = averageSensitivity;
compressedDev = deviationSensitivity;
compressedDev(isnan(compressedDev)) = [];
compressedAvg(isnan(compressedAvg)) = [];
RSE = 100*compressedDev./compressedAvg;
scatter(averageConcentration,RSE, 'k');
scatter(averageConcentration(RSE>15), RSE(RSE>15), 'r.');
grid on;
title('Residual Standard Error');
ylabel('RSE (%)');
xlabel('Concentration (mM)');

% Clarke Grid

logicalIndex = clarkeReal < 400;
clarkeReal = clarkeReal(logicalIndex);
clarkePredicted = clarkePredicted(logicalIndex);
logicalIndex = clarkePredicted < 400;
clarkePredicted = clarkePredicted(logicalIndex);
clarkeReal = clarkeReal(logicalIndex);
subplot(2,4,4);
[Count,Percent] = clarke(clarkeReal,clarkePredicted,900)
hold on;

% Recorded reference data from blood glucometer 
expected = [90 90 90 180 180 180 270 270 270 360 360 360 450 450 450,...
    90 90 90 150 150 150 210 210 210 270 270 270 330 330 330 390 390 390];
recorded = [93 99 95 205 204 207 307 300 313 378 382 388 446 430 440,...
    96 96 78 165 164 165 219 213 223 282 258 287 252 335 231 348 373 387];
scatter(expected,recorded,'ro','filled');

[~,~,RSQC,PC] = linearFit(clarkeReal,clarkePredicted)

[~,~,RSQG,PG] = linearFit(expected,recorded)

[RSQ1,P1] = oneOneLinearFit(clarkeReal,clarkePredicted)

[RSQ2,P2] = oneOneLinearFit(expected,recorded)

%% MARD Plot
output = figure(901);
set(output,'Position', [100, 100, 1000, 772]); % 8.5x11 ratio
colorKey1 = {'ko-','bo-','go-','co-','mo-','yo-','ro-'};

% Split sensitivity and concentration by hour/calibration run
%splitSensHour = splitVectorCat(summedSensitivity,summedHour);
%splitConcHour = splitVectorCat(summedConcentration,summedHour);
%splitElecHour = splitVectorCat(summedElectrodeNumber,summedHour);

% Split average values by calibration iteration and eliminate zero regression value
splitSensHour = splitVectorCat(calAverage,calIteration);
splitConcHour = splitVectorCat(calConcentration,calIteration);
if splitConcHour(1,1) == 0
    splitSensHour(1,:) = [];
    splitConcHour(1,:) = [];
end

% Perform regression on first hour/calibration run
if B == 1
    [~,~,~,P2] = quadraticFit(splitSensHour(:,1),splitConcHour(:,1));
else
    [~,~,~,P2] = linearFit(splitSensHour(:,1),splitConcHour(:,1));
end

n = 1;
% Cycle through and plot per-concentration MARD calculations for each calibration run
for i=1:1:size(splitSensHour,2)
    splitSensConc = splitVectorCat(splitSensHour(:,i),splitConcHour(:,i));
    concList = unique(splitConcHour(:,i)).';
    concList(isnan(concList)) = [];
    concExpected = polyval(P2,splitSensConc);
    concExpected(isnan(concExpected)) = [];
    
    for j=1:1:size(concExpected,2)
        for k=1:1:size(concExpected,1);
            mardAll(k,j) = 100*abs(concExpected(k,j)-concList(j))/abs(concList(j));
        end
    end
    
    for j=1:1:size(mardAll,2)
        mardAvg(j) = nanmean(mardAll(:,j));
        mardDev(j) = nanstd(mardAll(:,j));
    end
    
    figure(900);
    subplot(2,4,3);
    hold on; grid on;
    if i < 8
        errorbar(concList,mardAvg,mardDev,colorKey1{i});
    else
        errorbar(concList,mardAvg,mardDev,'o-','Color',colorListRandom(n,:));
    end
    figure(901);
    subplot(2,1,1);
    hold on; grid on;
    if i < 8
        errorbar(concList,mardAvg,mardDev,colorKey1{i});
    else
        errorbar(concList,mardAvg,mardDev,'o-','Color',colorListRandom(n,:));
        n = n+1; %last occurence of n in loop
    end
    overallMARD(i) = nanmean(mardAvg(isfinite(mardAvg)));
    
end

figure(900);
subplot(2,4,3);
title('MARD Per Concentration (Based on Initial Regression');
legend('Initial','Second','Third','Fourth','Fifth','Sixth','Seventh');
legend('Location','SouthwestOutside');
ylabel('MARD (%)');
xlabel('Concentration (mM)');

figure(901);
subplot(2,1,2);
title('MARD Per Concentration (Based on Initial Regression');
legend('Initial','Second','Third','Fourth','Fifth','Sixth','Seventh');
legend('Location','SouthwestOutside');
ylabel('MARD (%)');
xlabel('Concentration (mM)');
%x = 2:1:length(overallMARD);
%plot(x,overallMARD(2:length(overallMARD)),'ko--','LineWidth',2);
plot(overallMARD,'ko--','LineWidth',2);

title('Overall MARD By Calibration');
ylabel('MARD (%)');
xlabel('Calibration Number');

output = figure(901);
set(output,'PaperUnits','inches','PaperPosition',[0 0 12.96 8.64],'PaperSize',[12.96 8.64]);
if B == 1
    saveas(output,[outputFolder,testName,' - ',testName,' - MARD quadratic.pdf'],'pdf') %Save figure
else
    saveas(output,[outputFolder,testName,' - ',testName,' - MARD linear.pdf'],'pdf') %Save figure
end

% Upper left data panel - text information
figure(900);
ax = subplot(2,4,1);
text(0.5,1,['Calibration Plot for ', testAnalyte{1}, ' Sensors'],...
    'HorizontalAlignment', 'center', 'FontSize', 16);
text(0.5,0.85,testName,...
    'HorizontalAlignment', 'center', 'FontSize', 14 , 'Color', 'b');
text(0.5,0.65,'Initial Calibration Statistics',...
    'HorizontalAlignment', 'center', 'FontSize', 14 , 'Color', 'k');
if LODc > 0.5 %mM
    text(0.5,0.4,['LOD: ', num2str(LOD), ' nA'],...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
else
    text(0.5,0.4,['LOD: ', num2str(LOD), ' nA'],...
        'HorizontalAlignment', 'center', 'FontSize', 12);
end
if RSQ < 0.95 % R2 value
    text(0.5,0.3,['R2: ', num2str(RSQ)],...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
else
    text(0.5,0.3,['R2: ', num2str(RSQ)],...
        'HorizontalAlignment', 'center', 'FontSize', 12);
end
if sensitivitySlope < 0.1 % nA/mM
    text(0.5,0.2,['Sensitivity: ', num2str(sensitivitySlope), ' nA/mM'],...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
else
    text(0.5,0.2,['Sensitivity: ', num2str(sensitivitySlope), ' nA/mM'],...
        'HorizontalAlignment', 'center', 'FontSize', 12);
end
set(ax,'visible','off');

% Save figure 900 after adding MARD calculation
output = figure(900);
subplot(2,4,1);

%add MARD to overall plot (without using initial regression)
text(0.5,0.1,['Chip MARD: ', num2str(nanmean(overallMARD(2:length(overallMARD)))), ' %'],...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k');

set(output,'PaperUnits','inches','PaperPosition',[0 0 12.96 8.64],'PaperSize',[12.96 8.64]);
if B == 1
    saveas(output,[outputFolder,testName,' - ',testName,' quadratic.pdf'],'pdf') %Save figure
else
    saveas(output,[outputFolder,testName,' - ',testName,' linear.pdf'],'pdf') %Save figure
end

if W == 1
    output = figure(902);
    set(output,'PaperUnits','inches','PaperPosition',[0 0 12.96 8.64],'PaperSize',[12.96 8.64]);
    if B == 1
        saveas(output,[outputFolder,testName,' - ',testName,' quadratic even.pdf'],'pdf') %Save figure
    else
        saveas(output,[outputFolder,testName,' - ',testName,' linear even.pdf'],'pdf') %Save figure
    end
    
    % Upper left data panel - text information
    ax = subplot(2,4,1);
    text(0.5,1,['Calibration Plot for ', testAnalyte{1}, ' Sensors'],...
        'HorizontalAlignment', 'center', 'FontSize', 16);
    text(0.5,0.85,testName,...
        'HorizontalAlignment', 'center', 'FontSize', 14 , 'Color', 'b');
    text(0.5,0.65,'Initial Calibration Statistics',...
        'HorizontalAlignment', 'center', 'FontSize', 14 , 'Color', 'k');
    if LODc > 0.5 %mM
        text(0.5,0.4,['LOD: ', num2str(LOD), ' nA'],...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
    else
        text(0.5,0.4,['LOD: ', num2str(LOD), ' nA'],...
            'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    if RSQ < 0.95 % R2 value
        text(0.5,0.3,['R2: ', num2str(RSQ)],...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
    else
        text(0.5,0.3,['R2: ', num2str(RSQ)],...
            'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    if sensitivitySlope < 0.1 % nA/mM
        text(0.5,0.2,['Sensitivity: ', num2str(sensitivitySlope), ' nA/mM'],...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
    else
        text(0.5,0.2,['Sensitivity: ', num2str(sensitivitySlope), ' nA/mM'],...
            'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    set(ax,'visible','off');
end

end
%}
%% Global sensitivity plotting function
function [averageSensitivity,deviationSensitivity,WLfit,CLfit,RSQ,sensitivitySlope,...
    calAverage,calConcentration,calDeviation,calIteration,averageConcentration] = ...
    summaryPlot(splitConcElectrodeNumber,splitSensElectrodeNumber,splitHourElectrodeNumber,...
    colorList1,colorList2,colorListRandom)
    
global B V

    % Plot individual sensor data, calculate and plot 
for i=1:1:size(splitConcElectrodeNumber,2)
    if isnan(splitConcElectrodeNumber(1,i)) == 0
        %if size(splitConcElectrodeNumber,2) > 1
        splitConcentration = splitVectorCat(splitConcElectrodeNumber(:,i),splitHourElectrodeNumber(:,i));
        splitSensitivity = splitVectorCat(splitSensElectrodeNumber(:,i),splitHourElectrodeNumber(:,i));
        
        % Select for quadratic or linear regression
        if B == 1
            [WLfit,CLfit,~,~] = quadraticFit(splitConcentration(:,1),splitSensitivity(:,1));
        else
            [WLfit,CLfit,~,~] = linearFit(splitConcentration(:,1),splitSensitivity(:,1));
        end
        if CLfit(1) > CLfit(length(CLfit)) %Label lines red if inverted
            scatter(splitConcentration(:,1),splitSensitivity(:,1),'r');
            plot(WLfit,CLfit,'r');
        else
            scatter(splitConcentration(:,1),splitSensitivity(:,1),'k');
            plot(WLfit,CLfit,'k');
        end
    end
end
title('Individual Sensors - Initial Calibration');
xlabel('Concentration (mM)');
if V == 4
    ylabel('Average Change in Current (nA)');
else
    ylabel('Average Sensor Current (nA)');
end

subplot(2,4,[6 7]);
hold on; grid on;

% Reconstruct summation vectors from (potentially) split components
summedConcentration = reshape(splitConcElectrodeNumber,...
    size(splitConcElectrodeNumber,1)*size(splitConcElectrodeNumber,2),1);
summedSensitivity = reshape(splitSensElectrodeNumber,...
    size(splitSensElectrodeNumber,1)*size(splitSensElectrodeNumber,2),1);
summedHour = reshape(splitHourElectrodeNumber,...
    size(splitHourElectrodeNumber,1)*size(splitHourElectrodeNumber,2),1);

summedConcentration(isnan(summedConcentration)) = [];
summedSensitivity(isnan(summedSensitivity)) = [];
summedHour(isnan(summedHour)) = [];

% Split current and concentration values by test/time value
% All three data streams must be vertical vectors with no NaN values
splitConcentration = splitVectorCat(summedConcentration,summedHour);
splitSensitivity = splitVectorCat(summedSensitivity,summedHour);

calIteration = NaN;
calAverage = NaN;
calDeviation = NaN;
calConcentration = NaN;

n = 1;
% Cycle through calibration iterations
for i=1:1:size(splitSensitivity,2)
    
    splitSensitivityConc = splitVectorCat(splitSensitivity(:,i),splitConcentration(:,i));
    averageConcentration = unique(splitConcentration).';
    
    for j=1:1:size(splitSensitivityConc,2)
        averageSensitivity(j) = nanmean(splitSensitivityConc(:,j));
        deviationSensitivity(j) = nanstd(splitSensitivityConc(:,j));
    end
    
    averageSensitivity(isnan(averageSensitivity)) = [];
    deviationSensitivity(isnan(deviationSensitivity)) = [];
    averageConcentration(isnan(averageConcentration)) = [];
    
    % Knock out zero concentration value if present
    if averageConcentration(1) == 0
        averageSensitivity(1) = [];
        deviationSensitivity(1) = [];
        averageConcentration(1) = [];
    end
    % Plot quadratic or linear regression based on averages over electrode
    % set for each concentration. If regression is the first, perform fit
    % and save polynomial coefficients for MARD calculation.
    if B == 1
        [WLfit,CLfit,RSQ,P] = quadraticFit(averageConcentration,averageSensitivity);
        [~,~,~,P2] = quadraticFit(averageSensitivity,averageConcentration);
        sensitivitySlope = P(2);
        if i == 1
            initialReg = P2;
        end
    else
        [WLfit,CLfit,RSQ,P] = linearFit(averageConcentration,averageSensitivity);
        [~,~,~,P2] = linearFit(averageSensitivity,averageConcentration);
        sensitivitySlope = P(1);
        if i == 1
            initialReg = P2;
        end
    end
    if i < 8
        errorbar(averageConcentration,averageSensitivity,deviationSensitivity,colorList1{i});
        plot(WLfit,CLfit,colorList2{i},'LineWidth',2);
    else
        errorbar(averageConcentration,averageSensitivity,deviationSensitivity,'o',...
            'Color',colorListRandom(n,:));
        plot(WLfit,CLfit,'--','LineWidth',2,'Color',colorListRandom(n,:));
        n = n+1;
    end
    
    iteration(1:length(averageSensitivity)) = i;
    calIteration = vertcat(calIteration,iteration.');
    calAverage = vertcat(calAverage,averageSensitivity.');
    calDeviation = vertcat(calDeviation,deviationSensitivity.');
    calConcentration = vertcat(calConcentration,averageConcentration.');

end

% Remove any nan values from average calibrations (all electrodes pooled together)
calIteration(isnan(calIteration)) = [];
calAverage(isnan(calAverage)) = [];
calDeviation(isnan(calDeviation)) = [];
calConcentration(isnan(calConcentration)) = [];

% Group regression plot and legend information
grid on;
title('Group Regressions by Time');
labels = {'Initial Data','Initial Fit','Second Data','Second Fit','Third Data',...
    'Third Fit','Fourth Data','Fourth Fit','Fifth Data','Fifth Fit','Sixth Data',...
    'Sixth Fit','Seventh Data','Seventh Fit'};
h_legend = legend(labels,'Location','SouthEastOutside');
set(h_legend,'FontSize',8);
xlabel('Concentration (mM)');
ylabel('Sensor Current (nA)');

end

%% Sensitivity data separation
function [sensTime,sensData,sensConcentration,sensIteration,...
    sensInterval,sensLength,sensAnalyte,sensHour] = sensTestSplit(TimePoint,...
    Measurement,testConcentration,testIteration,testInterval,...
    testLength,testAnalyte,testType,testHour)

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

%% Sensitivity Plotting Function
function sensPlotting(graphIndex,sensorID,TimePoint,Measurement,...
    concentration,testLength,testInterval,testIteration,testAnalyte,testHour)

% Current scaling factor (1e9 = plot data in nanoamps, etc.)
scaleFactor = 1e9;

% Randomly generated color list for high hour iteration numbers
colorListRandom = [rand(1,3);rand(1,3);rand(1,3);rand(1,3);rand(1,3);rand(1,3);rand(1,3)];

% Global variables
global clarkeReal clarkePredicted sensitivitySlope
global B V
global summedSensitivity summedConcentration summedDeviation summedHour

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
    for j=1:1:length(hourListing)
        
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
        average(isnan(average)) = []
        deviation(isnan(deviation)) = []
        
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
                [WLfit,CLfit,RSQ(j),P] = quadraticFit(splitConcenSorted,average);
                [~,~,~,P2] = quadraticFit(average,splitConcenSorted);
                if j == 1
                    initialReg = P2;
                end
            else
                [WLfit,CLfit,RSQ(j),P] = linearFit(splitConcenSorted,average);
                [~,~,~,P2] = linearFit(average,splitConcenSorted);
                if j == 1
                    initialReg = P2;
                end
            end
        catch
            close all;
            error(['Issue reading Sensitivity text data. '...
                'Check Sensitivity file header and sizes.']);
        end
        
        % Extract coefficients for future use
        if B == 1
            sensitivitySlope(j) = P(2);
            sensitivityIntercept = P(3);
            clarkeSlope = P2(2);
            clarkeIntercept = P2(3);
            clarkeQuadratic = P2(1);
        else
            sensitivitySlope(j) = P(1);
            sensitivityIntercept = P(2);
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
            expectedEndAll = polyval(initialReg,splitspeed5*1e9);
        elseif V == 2
            expectedEndAll = polyval(initialReg,splitspeed10*1e9);
        elseif V == 3
            expectedEndAll = polyval(initialReg,splitspeed15*1e9);
        elseif V == 4
            expectedEndAll = polyval(initialReg,splitspeed30*1e9);
        elseif V == 5
            expectedEndAll = polyval(initialReg,splitspeedDev*1e9);
        else
            expectedEndAll = polyval(initialReg,splitspeedEnd*1e9);
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
    end
    
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
    LODc(j) = average(1)+deviation(1)*3;
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
    labels = {'Initial Data','Initial Fit','Second Data','Second Fit','Third Data',...
        'Third Fit','Fourth Data','Fourth Fit','Fifth Data','Fifth Fit','Sixth Data',...
        'Sixth Fit','Seventh Data','Seventh Fit'};
    h_legend = legend(labels,'Location','SouthEast');
    set(h_legend,'FontSize',8);
    
    title('Concentration vs. Current');
    if V == 5
        ylabel('Change in Measured Current (nA)');
    else
        ylabel('Measured Current (nA)');
    end
    xlabel([testAnalyte{1},' Concentration']);
    
    % Add text information to center panel of plot
    ax = subplot(2,4,1);
    text(0.5,1,['Calibration Plot for ', testAnalyte{1}, ' Sensor'],...
        'HorizontalAlignment', 'center', 'FontSize', 16);
    text(0.5,0.85,['Sensor ID: ', num2str(sensorID)],...
        'HorizontalAlignment', 'center', 'FontSize', 14 , 'Color', 'b');
    text(0.5,0.65,'Initial Calibration Statistics',...
        'HorizontalAlignment', 'center', 'FontSize', 14 , 'Color', 'k');
    if LODc > 0.1 %0.5 %mM
        text(0.5,0.4,['LOD: ', num2str(LODc), ' nA'],...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
    else
        text(0.5,0.4,['LOD: ', num2str(LODc), ' nA'],...
            'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    if RSQ < 0.95 % R2 value
        text(0.5,0.3,['R2: ', num2str(RSQ)],...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
    else
        text(0.5,0.3,['R2: ', num2str(RSQ)],...
            'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    if sensitivitySlope < 0.1 % nA/mM
        text(0.5,0.2,['Sensitivity: ', num2str(sensitivitySlope), ' nA/mM'],...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
    else
        text(0.5,0.2,['Sensitivity: ', num2str(sensitivitySlope), ' nA/mM'],...
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

end

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
    xticklabel_rotate(1:length(labels),90,labels);
    
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