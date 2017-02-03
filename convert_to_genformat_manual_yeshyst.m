function [MeasurementC_single,TimePointC_single,testIterationC_single,testConcentrationC_single,testIntervalC_single,testLengthC_single,testHourC_single,testTypeC_single,testAnalyteC_single, testInterferentC_single]...
    =convert_to_genformat_manual_yeshyst(sampled_time_data,sampled_current_data,gluc_incr,flush_num,max_concentration,calib_num,analyte,interferent,test_type)


%doesnt work with multiple flushes yet and needs to be doubleclicked
MeasurementC_single=sampled_current_data;

TimePointC_single=[];
testIterationC_single=[];
testConcentrationC_single=[];
testIntervalC_single=[];
testLengthC_single=[];
testTypeC_single={};
testAnalyteC_single={};
testInterferentC_single={};
calib_glucose_unique=[];
calib_glucose_single_calib=[0];
calib_glucose_merged=[];

%Build calibration vector hysterisis
uniq_points_per_calib = round(max_concentration/gluc_incr);

for i=1:uniq_points_per_calib
    calib_glucose_unique = [calib_glucose_unique; gluc_incr*i];
end
for i=1:flush_num
    calib_glucose_single_calib=[calib_glucose_single_calib; calib_glucose_unique];
    calib_glucose_single_calib=sort(calib_glucose_single_calib);    
end
direction = 1;
for i=1:ceil(size(MeasurementC_single,1)/length(calib_glucose_single_calib))
    if direction == 1
        calib_glucose_merged = [calib_glucose_merged; calib_glucose_single_calib];
    elseif direction == -1
        calib_glucose_merged = [calib_glucose_merged; flipud(calib_glucose_single_calib)];
    end
    direction = direction * -1;
end


%Build rest of data vectors
for i=1:size(sampled_current_data,1)
    TimePointC_single=vertcat(TimePointC_single,sampled_time_data(i,:)-sampled_time_data(i,1));
    testIterationC_single=vertcat(testIterationC_single,mod(i-1,flush_num)+1);
    testConcentrationC_single = vertcat(testConcentrationC_single,calib_glucose_merged(i));
    testIntervalC_single = vertcat(testIntervalC_single,sampled_time_data(1,2)-sampled_time_data(1,1));
    testLengthC_single = vertcat(testLengthC_single,TimePointC_single(end,end));
    testTypeC_single=[testTypeC_single; test_type];
    testAnalyteC_single=[testAnalyteC_single; analyte]; 
    testInterferentC_single=[testInterferentC_single; interferent];
    
end

%Build Calibration vector (each hysteresis cycle is 2 calibrations)
testHourC_single = zeros(size(MeasurementC_single,1),1);
for i=1:length(calib_glucose_single_calib):size(MeasurementC_single,1)
    testHourC_single(i:i+length(calib_glucose_single_calib)-1) = calib_num;
    calib_num = calib_num + 1;
end

%{
if hyst_yesorno == 1
    for i=calib_num:-1:1
        sampled_current_data_split = [sampled_current_data_split(1:length(calib_glucose_unique2)*i,:); sampled_current_data_split(length(calib_glucose_unique2)*i:end,:)];
        sampled_time_data_split = [sampled_time_data_split(1:length(calib_glucose_unique2)*i,:); sampled_time_data_split(length(calib_glucose_unique2)*i:end,:)];
        testinterval = [testinterval(1:length(calib_glucose_unique2)*i); testinterval(length(calib_glucose_unique2)*i:end)];
        testiter = [testiter(1:length(calib_glucose_unique2)*i); testiter(length(calib_glucose_unique2)*i:end)];
        calib_current = [calib_current(1:length(calib_glucose_unique2)*i,:); calib_current(length(calib_glucose_unique2)*i:end,:)];
    end
end
%}

end