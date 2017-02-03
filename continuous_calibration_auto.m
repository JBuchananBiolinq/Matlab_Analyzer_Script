function [sampled_time_no_offset,sampled_current,concentration,iteration,test_length,test_hour,test_type,test_analyte,test_interference,test_interval]=continuous_calibration_auto(time,current,flushnum,glucose_increment,max_glucose,time_between_increment,sample_length,graph_index,analyte,buffer_time,every_x_calibs, zero_offset_t, last_calib_num)

%variable initialization
concentration =[];
iteration=[];
test_length=[];
test_hour=[];
test_type={};
test_interval=[];
test_analyte={};
test_interference={};
sampled_time_no_offset=[];

%gather basic data info
points_per_calib=round(max_glucose/glucose_increment)+1;
total_points=round((time(end)+buffer_time)/(time_between_increment+buffer_time));
total_points=total_points-mod(total_points,points_per_calib);
time_resolution = time(2)-time(1);
num_calib=ceil(total_points/points_per_calib);
calibs_to_keep = 0:every_x_calibs:num_calib-1;

%Adjusts the number of calibrations that will be looked at based on user
%input (ie user can choose if they dont want to look at the last x number
%of calibrations).
if last_calib_num < num_calib
    num_calib = round(last_calib_num);
    total_points = num_calib * points_per_calib; 
end


%grab samples of current at precise intervals
for i =1:total_points   
    if any(floor((i-1)/points_per_calib)==calibs_to_keep)
        if i==1
            sampled_current = [current(round((time_between_increment-sample_length+zero_offset_t)/time_resolution):round((time_between_increment+zero_offset_t)/time_resolution))];
            sampled_time = [time(round((time_between_increment-sample_length+zero_offset_t)/time_resolution):round((time_between_increment+zero_offset_t)/time_resolution))];
        elseif round((time_between_increment*i+buffer_time*(i-1)+zero_offset_t)/time_resolution) <= length(time)     
            sampled_current = [sampled_current; current(round((time_between_increment*i+buffer_time*(i-1)-sample_length+zero_offset_t)/time_resolution):round((time_between_increment*i+buffer_time*(i-1)+zero_offset_t)/time_resolution))];
            sampled_time = [sampled_time; time(round((time_between_increment*i+buffer_time*(i-1)-sample_length+zero_offset_t)/time_resolution):round((time_between_increment*i+buffer_time*(i-1)+zero_offset_t)/time_resolution))];
        else
            sampled_current = [sampled_current; current(end-round(sample_length/time_resolution):end)];
            sampled_time = [sampled_time; time(end-round(sample_length/time_resolution):end)];
        end
        
        %stores data to be returned in the correct format for further
        %processing 
        concentration=[concentration; glucose_increment*floor(mod(i-1,points_per_calib)/flushnum)];
        iteration=[iteration; mod(i-1,flushnum)+1];
        test_interval=[test_interval;time_resolution];
        test_length=[test_length;sample_length];
        test_hour=[test_hour;floor((i-1)/(points_per_calib))];
        test_type=[test_type;'Sensitivity'];
        test_analyte=[test_analyte;analyte];
        test_interference=[test_interference; '[NaN]'];
        sampled_time_no_offset=[sampled_time_no_offset;sampled_time(end,:)-sampled_time(end,1)];
    end
end

%plot all data and then sampled data
figure(400+graph_index)
hold on
plot(time,current)
plot(sampled_time,sampled_current,'*g')
hold off

end