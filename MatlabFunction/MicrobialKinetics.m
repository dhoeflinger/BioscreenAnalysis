function [ lag_time, max_spec_growth_rate, OD_max, doubling_time, note ] = MicrobialKinetics(OD_values, time_interval, threshold, model)
%UNTITLED Summary of this function goes here
%  time interval assumed to be half hour time blocks



%set to zero initially
max_spec_growth_rate = 0;
note = '';
max_point = 0;
lag_time = 0;
doubling_time = 0;


%find the max OD value
[OD_max, index] = max(OD_values);
%find the min OD value 
[OD_min, min_index] = min(OD_values);



    
absolute_no_growth_threshold = threshold;

relative_no_growth_threshold = threshold + OD_min;

if (OD_max < absolute_no_growth_threshold && OD_max < relative_no_growth_threshold)
    note = 'No Growth Detected, Check Plot';    
    return;
end

if (OD_max < absolute_no_growth_threshold)
    note = 'Growth did not reach absolute threshold';
end

if (OD_max < relative_no_growth_threshold)
    note = 'Growth did not reach relative threshold';
end


non_log_slope =0;
[lag_time max_spec_growth_rate] = FindRegressionCurve(OD_values,time_interval, model);


%this is the old way
% for i=1:min(index+4, size(OD_values))
% 
%     if (i > 4) % if I have a time point 2 hours behind this one
%         if (OD_values(i-4) > OD_min && OD_values(i) > min(relative_no_growth_threshold, absolute_no_growth_threshold)) %and both OD reading is above 0 (to not screw up the log calc)
%           %find the slope at each time point with a two hour interval     
%           slope  =  (log(OD_values(i) - OD_min) - log(OD_values(i-4)- OD_min)) / (time_interval*4);
%        
%           
%           %if the slope is greater than the current max, set it to be the
%           %max and record the current time as the lag time
%           if (slope > max_spec_growth_rate)
%               max_spec_growth_rate = slope;
%               non_log_slope = (OD_values(i) - OD_values(i-4)) / (time_interval*4);
% 
%               max_point = (i-2);
%           end         
%         end
%     end
% end
% 
% if (max_point > 0)
%     lag_time = (max_point - ((OD_values(max_point) - OD_values(1)) / non_log_slope)) * time_interval;
% else
%     lag_time = 0;
% end




%use the max specific growth rate to calc the doubling time t_d = ln(2)/u
doubling_time = log(2) / max_spec_growth_rate;

lag_time_str = sprintf('lag time = %f, max growth = %f, doubling time = %f', lag_time, max_spec_growth_rate, doubling_time);

legend(lag_time_str, 'location', 'SouthOutside');
