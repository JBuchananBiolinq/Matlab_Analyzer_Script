function [MeasurementC_single,TimePointC_single,testIterationC_single,testConcentrationC_single,testIntervalC_single,testLengthC_single,testHourC_single,testTypeC_single,testAnalyteC_single, testInterferentC_single, StabilizeTime_increment_single]...
    =convert_to_genformat_manual_nohyst(sampled_time_data,sampled_current_data,test_type,gluc_incr,flush_num,calib_num,analyte,interferent,sel_incr,sel_interferents,max_gluc,stabilize_time)

MeasurementC_single=sampled_current_data;
increments_per_calib=max_gluc/gluc_incr+1;

TimePointC_single=[];
testIterationC_single=[];
testConcentrationC_single=[];
testIntervalC_single=[];
testLengthC_single=[];
testHourC_single=[];
testTypeC_single={};
testAnalyteC_single={};
testInterferentC_single={};


for i=1:size(MeasurementC_single,1)
    TimePointC_single = vertcat(TimePointC_single,sampled_time_data(i,:)-sampled_time_data(i,1));
    testIterationC_single = vertcat(testIterationC_single,mod(i-1,flush_num)+1);
    testIntervalC_single = vertcat(testIntervalC_single,sampled_time_data(1,2)-sampled_time_data(1,1));
    testLengthC_single = vertcat(testLengthC_single,TimePointC_single(end,end));
    testHourC_single = vertcat(testHourC_single,calib_num+floor(floor((i-1)/flush_num)/increments_per_calib));
    testTypeC_single = [testTypeC_single; test_type];
    testAnalyteC_single = [testAnalyteC_single;analyte];
       
    if strcmp(test_type,'Sensitivity')
        testInterferentC_single=[testInterferentC_single;interferent];
        testConcentrationC_single = vertcat(testConcentrationC_single, floor(mod((i-1),increments_per_calib)/(flush_num))*gluc_incr);  
    elseif strcmp(test_type,'Selectivity') && i == 1
        testInterferentC_single = [testInterferentC_single;'Baseline';'Control';sel_interferents];   
        testConcentrationC_single = [0;sel_incr;sel_incr];
    end
               
end

if strcmp(test_type,'Sensitivity')
    remove_stabtime_ind = find(testConcentrationC_single == 0)-1;
    remove_stabtime_ind(1) = [];
    stabilize_time(remove_stabtime_ind) = [];   
    StabilizeTime_increment_single = reshape(stabilize_time, increments_per_calib-1, length(stabilize_time)/(increments_per_calib-1))';
elseif strcmp(test_type,'Selectivity')
    StabilizeTime_increment_single = stabilize_time';
end

end


