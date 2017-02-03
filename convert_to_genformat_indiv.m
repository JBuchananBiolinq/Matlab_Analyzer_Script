function [MeasurementC_single,TimePointC_single,testIterationC_single,testConcentrationC_single,testIntervalC_single,testLengthC_single,testHourC_single,testTypeC_single,testAnalyteC_single, testInterferentC_single]=convert_to_genformat_manual_nohyst(sampled_time_data,sampled_current_data,test_type,gluc_incr,flush_num,calib_num,analyte,interferent)

MeasurementC_single=sampled_current_data;
TimePointC_single=sampled_time_data-sampled_time_data(1);

testIterationC_single=[];
testConcentrationC_single=[];
testIntervalC_single=[];
testLengthC_single=[];
testTypeC_single={};
testAnalyteC_single={};
testInterferentC_single={};

for i=1:size(MeasurementC_single,1)
    testIterationC_single=vertcat(testIterationC_single,mod(i-1,flush_num)+1);
    testConcentrationC_single = vertcat(testConcentrationC_single; floor((i-1)/(flush_num))*gluc_incr);
    testIntervalC_single = vertcat(testIntervalC_single,test_type);
    testLengthC_single = vertcat(test_LengthC_single,TimePointC_single(end,end));
    testHourC_single = vertcat(testHourC_single,calib_num);
    testTypeC_single=[testTypeC_single; test_type];
    testAnalyteC_single=[testAnalyteC_single;analyte];
    testInterferentC_single=[testInterferentC_single;interferent];
end 











