%% Global sensitivity plotting function
function [averageSensitivity,deviationSensitivity,WLfit,CLfit,RSQ,sensitivitySlope,...
    calAverage,calConcentration,calDeviation,calIteration,averageConcentration] = ...
    summaryPlot(splitConcElectrodeNumber,splitSensElectrodeNumber,splitHourElectrodeNumber,...
    colorList1,colorList2,colorListRandom)
    

if size(splitConcElectrodeNumber,1) > 1
    
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

% Plot group calibration/regression based on individual sensors
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

else
    averageSensitivity = NaN;
    deviationSensitivity = NaN;
    WLfit = NaN;
    CLfit = NaN;
    RSQ = NaN;
    sensitivitySlope = NaN;
    calAverage = NaN;
    calConcentration = NaN;
    calDeviation = NaN;
    calIteration = NaN;
    averageConcentration = NaN;
    
end

end
