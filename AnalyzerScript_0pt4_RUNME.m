function AnalyzerScript_0pt3()
% 0.3 New major version: Attempted major restructuring of data handling for
% faster and simpler operation
% Added continuous calibration data extraction on 6/22/16

% Required files for operation:
% dataImport_0pt1.m
% channelSeparation_0pt1.m
% clarke.m
% splitVectorCat.m

clear all; close all;

% Global variable initialization & definement
global clarkeReal clarkePredicted googleDrive 
global B V
global summedSensitivity summedConcentration summedCounter
global summedHour summedElectrodeNumber
global summedInterferenceAverage summedInterferenceLabel summedInterferenceHour
global summedSampleCurrentC  summedTestHour summedTestType summedElectrodelabel summedSampleCurrentCSensitivity
summedCounter = 1;
summedTestType = {};
sensitivitySlope_indiv=[];
LOD_indiv=[];
LOL_indiv=[];
LOQ_indiv=[];
R2_indiv=[];
polyfit_coef = [];

% Establish main file directory
googleDrive = strrep(pwd,'Matlab Script','');
%googleDrive = 'C:\Users\Public\Dropbox\Analyzer Results\';
%googleDrive = 'C:\Users\JP\Dropbox\Analyzer Results\';

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

% Extract list of files from directory
calibrationList = dir([rawCalibDir,'\*.txt']);

% Choose Test type and all the associated parameters
A = menu('Choose Test Type?','Continuous(Manual)','Individual','Continuous(Automated)');

if A == 1 || A == 3
    flush_num = inputdlg('Enter the Number of Flushes (i.e. iterations) per Glucose Increment:');
    %channel_ignore_string = inputdlg('Enter the electrode numbers you wish to ignore (seperate numbers by comma):');
    %channel_ignore_num = cellfun(@str2num,strsplit(channel_ignore_string{:},','));
    sample_time = inputdlg('Enter the desired Sample Time (s):');
    glucose_increment = inputdlg('Enter the Glucose Increment (mM):');
    max_gluc_conc = inputdlg('Enter Maximum Glucose Concentration(mM):');
  
    if A == 1
        end_or_middle = menu('Choose Sample location','End','Middle');
        hyst_yesorno = menu('Will Hysteresis be Tested?','Yes','No');
        selectivity_yesorno = menu('Will selectivity data be Included?','Yes','No');
        same_boundary_ind = menu('Do glucose increments occur at the same time for all channels?','Yes','No');
    
        if selectivity_yesorno == 1
            selectivity_increment = inputdlg('Enter Glucose Increment for Selectivty Data(mM):');
            selectivity_interferents = inputdlg('Enter Interferents(ex: AA+UA+ACT):'); 
        else
             selectivity_increment = {'0'};
             selectivity_interferents = {'[NaN]'};
        end
     
    elseif A == 3
        time_between_increment = inputdlg('Enter the time between increments (s):');
        zero_offset_t = inputdlg('Enter the zero offset time (s):');
        last_calib_num = inputdlg('Enter the first X number of Calibrations you wish to look at:');
        buffer_time = inputdlg('Enter the buffer/pump time (s):');
        every_x_calibrations = inputdlg('Record every X calibrations:');
    end
end 

B = menu('Quadratic or Linear Fit?','Quadratic','Linear');
if A == 1
    W = 2;
elseif A == 2
    V = menu('Select Sample Time','5s','10s','15s','30s','Derivative Last 2s','60s/End');
    W = menu('Split odd and even channels?','Yes','No');
end
clarkeReal = NaN;
clarkePredicted = NaN;

% Import all data from files and names of files within selected folder
[TimePoint,~,Measurement,...
    testConcentration,testIteration,testInterval,testLength,...
    channelID,testType,testAnalyte,testInterferent,testHour]...
    = dataImport_0pt1(rawCalibDir,calibrationList);

%% Split channel data & individual sensor plotting - loop through each available electrode
channelIDlist = unique(channelID);

%complete data variable intialization
boundary_index=cell(1,length(testHour));
StabilizeTime_increment={};
calib_we_label={};

%Cycle through each of the electrodes
for i=1:1:length(channelIDlist)
    
    % Extract information for individual channel. Only Indivudal
    % calibrations are processed completely with this function
    [MeasurementC,TimePointC,testConcentrationC,testIterationC,...
        testIntervalC,testLengthC,testHourC,testTypeC,testAnalyteC,testInterferentC]...
        = channelSeparation_0pt1(channelIDlist(i),channelID,Measurement,TimePoint,testConcentration,...
        testIteration,testInterval,testLength,testHour,testType,testAnalyte,testInterferent); 
    
    % Perform analysis of continuous calibration files, plot individual
    % curves. This basically converts continuous data into a general format
    % that can be processed by all of the plotting functions
    
    if A == 1 || A == 3 
        %create copy of initial extracted data
        MeasurementC_copy=MeasurementC;
        TimePointC_copy= TimePointC;
        testHourC_copy= testHourC;
        testTypeC_copy= testTypeC;
        testAnalyteC_copy= testAnalyteC;
        testInterferentC_copy= testInterferentC;
    
        %reinitialize data variables
        MeasurementC=[];
        TimePointC=[];
        testConcentrationC=[];
        testIterationC=[];
        testIntervalC=[];
        testLengthC=[];
        testHourC=[];
        testTypeC={};
        testAnalyteC={};
        testInterferentC={};
    end
    
    if A == 1 % continuous manual file preprocessing to general format
        for j=1:size(MeasurementC_copy,1) %cycle through calibrations
                
                %extract basic data for a single calibration
                non_zero_ind=find(TimePointC_copy(j,:) ~= 0);
                [sampled_current_data, sampled_time_data,stabilize_time,boundary_index{j}] = ...
                    continuous_calibration_interactive(TimePointC_copy(j,non_zero_ind)', MeasurementC_copy(j,non_zero_ind)', str2num(flush_num{:}) ,str2num(glucose_increment{:}), end_or_middle,str2num(sample_time{:}),903+((length(channelIDlist))+2*size(MeasurementC_copy,1))*i+2*j,channelIDlist(i),testTypeC_copy(j),testHourC_copy(j),str2num(max_gluc_conc{:}),hyst_yesorno,boundary_index{j});               
                
                if same_boundary_ind == 2
                    boundary_index{j}=[]; %removes saved increment locations if each channel occurs under differnet conditions
                end
                
                if hyst_yesorno == 2 | strcmp(testTypeC_copy(j),'Selectivity')
                    [MeasurementC_single,TimePointC_single,testIterationC_single,testConcentrationC_single,testIntervalC_single,testLengthC_single,testHourC_single,testTypeC_single,testAnalyteC_single, testInterferentC_single,StabilizeTime_increment_single]...
                        =convert_to_genformat_manual_nohyst(sampled_time_data,sampled_current_data,testTypeC_copy(j),str2num(glucose_increment{:}),str2num(flush_num{:}),testHourC_copy(j),testAnalyteC_copy(j),testInterferentC_copy(j),str2num(selectivity_increment{:}),selectivity_interferents{:},str2num(max_gluc_conc{:}),stabilize_time);
                elseif hyst_yesorno == 1
                    [MeasurementC_single,TimePointC_single,testIterationC_single,testConcentrationC_single,testIntervalC_single,testLengthC_single,testHourC_single,testTypeC_single,testAnalyteC_single, testInterferentC_single]...
                        =convert_to_genformat_manual_yeshyst(sampled_time_data,sampled_current_data,str2num(glucose_increment{:}),str2num(flush_num{:}),str2num(max_gluc_conc{:}),testHourC_copy(j),testAnalyteC_copy(j),testInterferentC_copy(j),testTypeC_copy(j));
                end
                
                MeasurementC=vertcat(MeasurementC,MeasurementC_single);
                TimePointC=vertcat(TimePointC,TimePointC_single);
                testConcentrationC=vertcat(testConcentrationC,testConcentrationC_single);
                testIterationC=vertcat(testIterationC,testIterationC_single);
                testIntervalC=vertcat(testIntervalC,testIntervalC_single);
                testLengthC=vertcat(testLengthC,testLengthC_single);
                testHourC=vertcat(testHourC,testHourC_single);
                testTypeC=[testTypeC;testTypeC_single];
                testAnalyteC=[testAnalyteC;testAnalyteC_single];
                testInterferentC=[testInterferentC;testInterferentC_single];  
                summedSampleCurrentC = vertcat(summedSampleCurrentC,sampled_current_data);
                if isempty(strmatch('Selectivity',testTypeC_single))
                    summedSampleCurrentCSensitivity = vertcat(summedSampleCurrentCSensitivity,sampled_current_data);
                end
                
                for k=1:size(StabilizeTime_increment_single,1)
                    StabilizeTime_increment=[StabilizeTime_increment; StabilizeTime_increment_single(k,:)];
                end
                        
        end
        
    elseif A == 3 % continuous automated file preprocessing to general format
        non_zero_ind=find(TimePointC_copy(1,:) ~= 0);
        [TimePointC,MeasurementC,testConcentrationC,testIterationC,testLengthC,testHourC,testTypeC,testAnalyteC,testInterferentC,testIntervalC]=...
            continuous_calibration_auto(TimePointC_copy(1,non_zero_ind),MeasurementC_copy(1,non_zero_ind),str2num(flush_num{:}),str2num(glucose_increment{:}),str2num(max_gluc_conc{:}),str2num(time_between_increment{:}),str2num(sample_time{:}),i,testTypeC_copy(1),str2num(buffer_time{:}),str2num(every_x_calibrations{:}),str2num(zero_offset_t{:}));        
        %no need for further processing since all the calibrations are
        %present in a single file 
    end 
    
    %Build Global Variables used for data writing Labels
    summedTestHour = vertcat(summedTestHour,testHourC);
    summedTestType = vertcat(summedTestType,testTypeC);
    summedElectrodelabel = vertcat(summedElectrodelabel,i*ones(length(testHourC),1)); 
        
    % Perform full suite of extraction. Continuous files must be
    % preprocessed first 
    
        % Extract sensitivity data for individual channel
        [sensTime,sensData,sensConcentration,sensIteration,...
            sensInterval,sensLength,sensAnalyte,sensHour] = sensTestSplit(TimePointC,...
            MeasurementC,testConcentrationC,testIterationC,testIntervalC,...
            testLengthC,testAnalyteC,testTypeC,testHourC);
        
        % Plot sensitivity data for individual channel     
        
        [sensitivitySlope_temp,LOD_temp,LOL_temp,LOQ_temp,R2_temp,polyfit_coef_temp] = sensPlotting(i*10,channelIDlist(i),length(channelIDlist),sensTime,sensData,sensConcentration,...
            sensLength,sensInterval,sensIteration,sensAnalyte,sensHour,1+str2num(max_gluc_conc{:})/str2num(glucose_increment{:})); 
        sensitivitySlope_indiv = [sensitivitySlope_indiv; sensitivitySlope_temp];
        LOD_indiv = [LOD_indiv; LOD_temp];
        LOL_indiv = [LOL_indiv; LOL_temp];
        LOQ_indiv = [LOQ_indiv; LOQ_temp];
        R2_indiv = [R2_indiv; R2_temp];
        polyfit_coef = [polyfit_coef; polyfit_coef_temp];
        
        % Create/append electrode number vector for pooled data analysis
        if isnan(sensConcentration) == 0
            testHourC(isnan(testHourC))=[];
            electrodeNumber(1:length(sensConcentration)/length(unique(sensIteration))) = channelIDlist(i);
            summedElectrodeNumber = vertcat(summedElectrodeNumber,electrodeNumber.');
        end
        
        % Extract selectivity data for individual channel
        [selectTime,selectData,selectConcentration,selectIteration,...
            selectInterval,selectLength,selectAnalyte,selectInterferent,selectHour]...
            = selTestSplit(TimePointC,MeasurementC,testConcentrationC,testIterationC,...
            testIntervalC,testLengthC,testAnalyteC,testTypeC,testInterferentC,testHourC);
        
        % Plot selectivity data for individual channel
        selPlotting(i*10,channelIDlist(i),selectTime,selectData,...
            selectConcentration,selectLength,selectInterval,selectIteration,...
            selectAnalyte,selectInterferent,selectHour);
        
        % Extract stability data for individual channel
        [~,stabData,~,stabIteration,...
            stabInterval,stabLength,~,stabHour] =...
            stabTestSplit(TimePointC,MeasurementC,testConcentrationC,testIterationC,...
            testIntervalC,testLengthC,testAnalyteC,testTypeC,testHourC);
        
        % Plot stability data for individual channel
        stabPlotting(i*10,channelIDlist(i),stabData,stabHour,...
            stabLength,stabInterval,stabIteration);
        
        output = figure(i*10);
        set(output,'PaperUnits','inches','PaperPosition',[0 0 12.96 8.64],'PaperSize',[12.96 8.64]);
        saveas(output,[outputFolder,testName,' - ',num2str(channelIDlist(i)),'.pdf'],'pdf') %Save figure
end

%individual electrode calibrations summary graph save
output = figure(123);
set(output,'PaperUnits','inches','PaperPosition',[0 0 12.96 8.64],'PaperSize',[12.96 8.64]);
saveas(output,[outputFolder,testName,'_Individual_Calibrations_Summary.pdf'],'pdf') %Save figure

%Write Extracted Data to text file 
calib_we_label={}; %build table data labels

for i=1:length(summedElectrodelabel)
    calib_we_label{i}=strcat('WE',num2str(summedElectrodelabel(i)),'_Calib',num2str(summedTestHour(i)),'_',summedTestType{i},':');
end 

calib_we_label = unique(calib_we_label,'stable');

%write stabilize time data to text file
stabilize_time_table = table(calib_we_label',StabilizeTime_increment);
writetable(stabilize_time_table,strcat(outputFolder,testName,'_StabilizeTime.txt'),'Delimiter','\t');

%write calibration points to text file
calib_write_cell={};
calib_current_points = mean(summedSampleCurrentC ,2);

for i=1:length(calib_current_points)
    if i==1
        hour=summedTestHour(1);
        type=summedTestType{1};
        electrode_num=summedElectrodelabel(1);
        calib_write_cell{1}=[calib_current_points(1)];        
    elseif ((hour == summedTestHour(i) || (isnan(hour) == 1 && isnan(summedTestHour(i))==1))&& electrode_num == summedElectrodelabel(i) && strcmp(type,summedTestType{i}))|| i == length(calib_current_points)
        calib_write_cell{end}=[calib_write_cell{end}, calib_current_points(i)]; 
    else
        hour=summedTestHour(i);
        type=summedTestType{i};
        electrode_num=summedElectrodelabel(i);
        calib_write_cell{end+1}=[calib_current_points(i)]; 
    end
end

calib_points_table = table(calib_we_label',calib_write_cell');

writetable(calib_points_table,strcat(outputFolder,testName,'_CalibPoints.txt'),'Delimiter','\t');

%write sensitivity,LOD,LOL,LOQ,R2 values to text file
sense_labels ={};
LOx_labels ={};
R2_labels ={};
for i=1:length(unique(summedElectrodelabel))
    sense_labels{end+1} = strcat('Sensor ',num2str(i),' (nA/mM)');
    LOx_labels{end+1} = strcat('Sensor ',num2str(i),' (mM)');
    R2_labels{end+1} = strcat('Sensor ',num2str(i));
end

sens_table = table(sense_labels',num2cell(sensitivitySlope_indiv));
writetable(sens_table,strcat(outputFolder,testName,'_Sensitivity.txt'),'Delimiter','\t');

LOL_table = table(LOx_labels',num2cell(LOL_indiv));
writetable(LOL_table,strcat(outputFolder,testName,'_LOL.txt'),'Delimiter','\t');

LOQ_table = table(LOx_labels',num2cell(LOQ_indiv));
writetable(LOQ_table,strcat(outputFolder,testName,'_LOQ.txt'),'Delimiter','\t');

LOD_table = table(LOx_labels',num2cell(LOD_indiv));
writetable(LOD_table,strcat(outputFolder,testName,'_LOD.txt'),'Delimiter','\t');

R2_table = table(R2_labels',num2cell(R2_indiv));
writetable(R2_table,strcat(outputFolder,testName,'_R2.txt'),'Delimiter','\t');

%Write polyfit coefficient data to text file 
polyfit_coeff_table = table(R2_labels',polyfit_coef);
writetable(polyfit_coeff_table,strcat(outputFolder,testName,'_polyfit_coefficient.txt'),'Delimiter','\t');

if A == 2 || A == 1
%% Plot data for full set of electrodes (ie summary graphs)
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
    [averageSensitivity,deviationSensitivity,WLfit,CLfit,RSQ,TotalsensitivitySlope,...
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
    
    splitConcElectrodeNumber;
    splitSensElectrodeNumber;
    splitHourElectrodeNumber;
    
    [averageSensitivity,deviationSensitivity,WLfit,CLfit,RSQ,TotalsensitivitySlope,...
        calAverage,calConcentration,calDeviation,calIteration,averageConcentration] = summaryPlot(...
        splitConcElectrodeNumber,splitSensElectrodeNumber,...
        splitHourElectrodeNumber,colorList1,colorList2,colorListRandom);
end

LOD =(averageSensitivity(1) + deviationSensitivity(1)*3);
LOQ =(averageSensitivity(1) + deviationSensitivity(1)*10);
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
    totalLabel = {'',totalLabel{:},''}
    set(gca,'xticklabel',totalLabel)
    set(gca, 'FontSize', 10)
    %xticklabel_rotate(1:length(totalLabel),90,totalLabel);
else
    text(0.5,0.5,'No Selectivity Data',...
        'HorizontalAlignment', 'center', 'FontSize', 14);
    set(ax,'visible','off');
end

% RSD Plot
if size(averageConcentration,1) > 1
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
end


% Clarke Grid

if size(clarkeReal,1) > 1
    
logicalIndex = clarkeReal < 400;
clarkeReal = clarkeReal(logicalIndex);
clarkePredicted = clarkePredicted(logicalIndex);
logicalIndex = clarkePredicted < 400;
clarkePredicted = clarkePredicted(logicalIndex);
clarkeReal = clarkeReal(logicalIndex);
subplot(2,4,4);
[Count,Percent] = clarke(clarkeReal,clarkePredicted,900);
hold on;

% Recorded reference data from blood glucometer 
expected = [90 90 90 180 180 180 270 270 270 360 360 360 450 450 450,...
    90 90 90 150 150 150 210 210 210 270 270 270 330 330 330 390 390 390];
recorded = [93 99 95 205 204 207 307 300 313 378 382 388 446 430 440,...
    96 96 78 165 164 165 219 213 223 282 258 287 252 335 231 348 373 387];
scatter(expected,recorded,'ro','filled');

[~,~,RSQC,PC] = linearFit(clarkeReal,clarkePredicted);

[~,~,RSQG,PG] = linearFit(expected,recorded);

[RSQ1,P1] = oneOneLinearFit(clarkeReal,clarkePredicted);

[RSQ2,P2] = oneOneLinearFit(expected,recorded);

end
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
clc
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
    text(0.5,0.5,['LOD: ', num2str(LOD), ' nA'],...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
    text(0.5,0.4,['LOQ: ', num2str(LOQ), ' nA'],...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
else
    text(0.5,0.5,['LOD: ', num2str(LOD), ' nA'],...
        'HorizontalAlignment', 'center', 'FontSize', 12);
    text(0.5,0.4,['LOQ: ', num2str(LOQ), ' nA'],...
        'HorizontalAlignment', 'center', 'FontSize', 12);
end
if RSQ < 0.95 % R2 value
    text(0.5,0.3,['R2: ', num2str(RSQ)],...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
else
    text(0.5,0.3,['R2: ', num2str(RSQ)],...
        'HorizontalAlignment', 'center', 'FontSize', 12);
end
if TotalsensitivitySlope < 0.1 % nA/mM
    text(0.5,0.2,['Sensitivity: ', num2str(TotalsensitivitySlope), ' nA/mM'],...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
else
    text(0.5,0.2,['Sensitivity: ', num2str(TotalsensitivitySlope), ' nA/mM'],...
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
    if TotalsensitivitySlope < 0.1 % nA/mM
        text(0.5,0.2,['Sensitivity: ', num2str(TotalsensitivitySlope), ' nA/mM'],...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
    else
        text(0.5,0.2,['Sensitivity: ', num2str(TotalsensitivitySlope), ' nA/mM'],...
            'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    set(ax,'visible','off');
end
end
end