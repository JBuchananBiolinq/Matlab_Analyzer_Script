%% Continuous sensitivity plotting function
function sensContinuous(graphIndex,sensorID,average,deviation,concentration,testAnalyte,calibration_num,testHour)

% Global variables
global clarkeReal clarkePredicted sensitivitySlope
global B V
global summedSensitivity summedConcentration summedDeviation summedHour

scaleFactor = 1;

hourListing = unique(testHour);

for j=1:length(hourListing)

% Sort concentrations & values if out of order, remove NaN values
[splitConcenSorted,sortIndex] = sort(concentration(:,j));
%average = scaleFactor*average(sortIndex,j);
%deviation = scaleFactor*deviation(sortIndex,j);
splitConcenSorted(isnan(splitConcenSorted)) = [];
%average(isnan(average)) = [];
%deviation(isnan(deviation)) = [];

% Plot sensitivity measurements by calibration #
% elseif statements included for backwards compatibility
fig2 = figure(graphIndex);
set(fig2,'Position', [100, 200, 1800, 740]); % 8.5x11 ratio
subplot(2,4,[5 6]);
hold on; grid on;

colorList1 = {'ko','bo','go','co','mo','yo','ro'}; %errorbar
colorList2 = {'k--','b--','g--','c--','m--','y--','r--'}; %regression lines
colorList3 = {'k','b','g','c','m','y','r'}; %scatter plot
colorList4 = {'ko-','bo-','go-','co-','mo-','yo-','ro-'}; %errorbar

errorbar(splitConcenSorted,average(:,j),deviation(:,j),colorList1{j});

% Add calibration/sensitivity data to full data set vectors
hourList(:,1:length(average)) = hourListing(j);
if size(average,1) > 1
    summedSensitivity = vertcat(summedSensitivity,average(:,j));
else
    summedSensitivity = vertcat(summedSensitivity,average(:,j).');
end
summedDeviation = vertcat(summedDeviation,deviation(:,j));
summedConcentration = vertcat(summedConcentration,splitConcenSorted);
summedHour = vertcat(summedHour,hourList.');

% Calculate regression of sensitivity information (linear or quad)
%try
    if B == 1
        [WLfit,CLfit,RSQ(j),P] = quadraticFit(splitConcenSorted,average(:,j));
        [~,~,~,P2] = quadraticFit(average(:,j),splitConcenSorted);
        if j == 1
            initialReg = P2;
        end
    else
        [WLfit,CLfit,RSQ(j),P] = linearFit(splitConcenSorted,average(:,j));
        [~,~,~,P2] = linearFit(average(:,j),splitConcenSorted);
        if j == 1
            initialReg = P2;
        end
    end
%catch
%    close all;
%    error(['Issue reading Sensitivity text data. '...
%        'Check Sensitivity file header and sizes.']);
%end

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
plot(WLfit,CLfit,colorList2{j},'LineWidth',2);


% Calculate and plot residual standard error values
% elseif statements for backwards compatibility
subplot(2,4,2);
hold on;
concentrationList = unique(concentration(:,j));
RSE = 100*deviation(:,j)./average(:,j);
scatter(concentrationList,RSE,colorList3{j});


% If values over 5%, mark with solid red circle
scatter(concentrationList(RSE>5), RSE(RSE>5), 'r.');
grid on;
title('Residual Standard Error');
ylabel('RSE (%)');
xlabel('Concentration (mM)');

% MARD Plot - calculate and plot by concentration
subplot(2,4,3);
hold on; grid on;

expectedEndAll = polyval(initialReg,average);

for i=1:1:length(concentrationList)
    mardAll(i) = 100*abs(nanmean(expectedEndAll(i))-concentrationList(i))/abs(concentrationList(i));
end

%if hourListing(j) < 7
    plot(concentrationList,mardAll,colorList4{j});
%else
    %plot(concentrationList,mardAll(j,:),'Color',colorListRandom(n,:));
    %n = n+1; % Last graph, increase random color index by 1
%end

title('MARD Per Concentration');
ylabel('MARD (%)');
xlabel('Concentration (mM)');

subplot(2,4,4);
hold on; grid on;

for j=1:1:size(mardAll,1)
    mardAvg(j) = nanmean(mardAll(j,2:size(mardAll,2)));
    mardDev(j) = nanstd(mardAll(j,2:size(mardAll,2)));
end

errorbar(1,mardAvg,mardDev,'k--','LineWidth',2);

title('Average MARD Per Calibration');
ylabel('MARD (%)');
xlabel('Calibration Number');

% Aggregate and approximate Clarke error grid values from linear regression
aggconcen = NaN;
aggspeed = NaN;
%for i=1:1:size(splitconcen,2)
%    aggconcen = vertcat(aggconcen,splitconcen(:,i));
%    aggspeed = vertcat(aggspeed,splitspeedEnd(:,i));
%end

%clarkeReal = vertcat(clarkeReal,18*aggconcen);
%clarkePredicted = vertcat(clarkePredicted,18.*((clarkeQuadratic*(1e9*aggspeed).^2)...
%    +(clarkeSlope*1e9*aggspeed)+clarkeIntercept));

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

end

% Add legend and axis data to calibration/sensititivity plot
subplot(2,4,[5 6]);
labels = {'Initial Data','Initial Fit','Second Data','Second Fit','Third Data',...
    'Third Fit','Fourth Data','Fourth Fit','Fifth Data','Fifth Fit','Sixth Data',...
    'Sixth Fit','Seventh Data','Seventh Fit'};
h_legend = legend(labels,'Location','SouthEast');
set(h_legend,'FontSize',8);

title('Concentration vs. Current');
ylabel('Measured Current (nA)');
xlabel([testAnalyte{1},' Concentration']);

% Add text information to center panel of plot
ax = subplot(2,4,1);
text(0.5,1,['Calibration Plot for ', testAnalyte{1}, ' Sensor'],...
    'HorizontalAlignment', 'center', 'FontSize', 16);
text(0.5,0.85,['Sensor ID: ', num2str(sensorID)],...
    'HorizontalAlignment', 'center', 'FontSize', 14 , 'Color', 'b');
text(0.5,0.65,['Calibration ' num2str(calibration_num) -1 ' Statistics'],...
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

end