function [sampled_current_data_split,sampled_time_data_split,stabilize_time,boundary_index]=continuous_calibration_interactive(time_data,current_data,flush_num,gluc_incr,test_measure_loc,measure_interval,figure_num,channelnum,test_type,calib_num,max_concentration,hyst_yesorno,presaved_boundary_index)

%Variable Initialization
measure_range=[];
stabilize_time=[];
sampled_time_data=[];
sampled_current_data=[];
stabilize_indices=[];
%calib_std=[0];
time_increment=time_data(2)-time_data(1);
time_data=time_data(find(~isnan(current_data)));
current_data=current_data(find(~isnan(current_data)));


%Graph raw data
fh_fig = figure(figure_num);
hold on 
plot(time_data(1:end-1)/60,current_data(1:end-1),'b');
title([test_type ' Calibration of Raw Data Points for channel ' num2str(channelnum)]);
xlabel('Time (min)');
ylabel('Current (A)');

NumTicks = 30;
set(gca,'FontSize',7)
L = get(gca,'XLim');
set(gca,'XMinorTick','on');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

    
if isempty(presaved_boundary_index) 
    %grab initial boundary estimates from user
    [boundary_estimates_x,boundary_estimates_y] = getpts(fh_fig);
    boundary_estimates_x = boundary_estimates_x * 60;
    boundary_index_init=dsearchn([time_data/max(time_data),current_data/max(current_data)],[boundary_estimates_x/max(time_data),boundary_estimates_y/max(current_data)]);
    boundary_index_init=sort(boundary_index_init);
    num_calibs_perfile = length(boundary_index_init)/(max_concentration/gluc_incr);

    %filtering for final boundary points
    boundary_index=[boundary_index_init(1)];
    for i=2:length(boundary_index_init)
    
        if boundary_index_init(i)+round(time_data(end)*0.01/num_calibs_perfile) < length(current_data)
            test_range=[boundary_index_init(i)-round(time_data(end)*0.01/num_calibs_perfile):boundary_index_init(i)+round(time_data(end)*0.01/num_calibs_perfile)];
        else
            test_range=[boundary_index_init(i)-round(time_data(end)*0.01/num_calibs_perfile):length(current_data)];
        end
        
        test_range = test_range(find(test_range > 0)); %test indices must be greater than 0
        current_changes=abs(diff(current_data(test_range))/mean(current_data(test_range)));   
    
        if isempty(find(current_changes > 0.02,1))
            boundary_index(i)=test_range(end);
        else
            boundary_index(i)=find(current_changes > 0.02,1)+test_range(1)-1;
            if abs(current_data(boundary_index(i))-current_data(boundary_index_init(i)))/mean(current_data)>0.05
                boundary_index(i)=boundary_index_init(i);
            end 
        end
    end
    
else    
    boundary_index = presaved_boundary_index;
end

%find maximum end currents at the boundary of each flushing
for i=1:length(boundary_index)-1
    [~,max_current_end_index(i)] = max(current_data(boundary_index(i+1)-50:boundary_index(i+1)));
    max_current_end_index(i) =  max_current_end_index(i) + boundary_index(i+1)-51;
    max_current_end_current(i) = mean(current_data(max_current_end_index(i)-20:max_current_end_index(i)));
    current_data(max_current_end_index(i)-10:max_current_end_index(i));
    mean(current_data(max_current_end_index(i)-10:max_current_end_index(i)));
end 

%Zero glucose calibration point
zero_range=[round(boundary_index(1)/2):round(boundary_index(1)/2)+measure_interval/(time_increment)];
calib_current=[mean(current_data(zero_range))];
sampled_current_data_split = current_data(zero_range)';
sampled_time_data_split = time_data(zero_range)';

%find all non-zero calibration points
for j=0:(length(boundary_index)-1)/flush_num-1
    for i=1:flush_num
        if test_measure_loc == 2
            measure_range=[measure_range, round(median(boundary_index(flush_num*j+i):1:boundary_index(flush_num*j+i+1))-0.5*measure_interval/time_increment),round(median(boundary_index(flush_num*j+i):1:boundary_index(flush_num*j+i+1))+0.5*measure_interval/time_increment)];
        elseif test_measure_loc ==1
            measure_range=[measure_range, boundary_index(flush_num*j+i+1)-measure_interval/time_increment,boundary_index(flush_num*j+i+1)];
        end 
        
        %Find stabilization points
        for k=boundary_index(flush_num*j+i)+30:1:max_current_end_index(flush_num*j+i)
            if current_data(k)- calib_current(1)> 0.9*(max_current_end_current(flush_num*j+i)-calib_current(1)) && abs(mean(current_data(k-50:k+50))-current_data(k))/mean(current_data(k-30)) < 0.05
                stabilize_indices = [stabilize_indices k];
                break
            elseif k == max_current_end_index(flush_num*j+i)
                stabilize_indices = [stabilize_indices k];
            end
        end
        
        if length(stabilize_indices) > 0
            stabilize_time=[stabilize_time; time_increment*(stabilize_indices(end)-boundary_index(flush_num*j+i))];
        end
        
        sample_avg(i)=mean(current_data(measure_range(end-1):measure_range(end)));
        
        sampled_time_data=vertcat(sampled_time_data,time_data(measure_range(end-1):measure_range(end)));
        sampled_time_data_split=vertcat(sampled_time_data_split,[time_data(measure_range(end-1):measure_range(end))]');
        sampled_current_data=vertcat(sampled_current_data, current_data(measure_range(end-1):measure_range(end)));
        sampled_current_data_split = vertcat(sampled_current_data_split,[current_data(measure_range(end-1):measure_range(end))]');
    end 
    %calib_current=[calib_current; mean(sample_avg)];
    %calib_std=[calib_std; std(sample_avg)];
end

%plotting
plot(vertcat(time_data(zero_range),sampled_time_data)/60,vertcat(current_data(zero_range),sampled_current_data),'*g');
plot(time_data(boundary_index)/60,current_data(boundary_index),'ro')
plot(time_data(stabilize_indices)/60,current_data(stabilize_indices),'ko')
legend('Raw Data','Sampled Data','Boundaries','Stabilization Point','Location','northwest');

for i=1:length(stabilize_indices)
    text(time_data(stabilize_indices(i))/60,min(current_data),strcat(num2str(stabilize_time(i)),' s'),'FontSize',6);
end

hold off
end


