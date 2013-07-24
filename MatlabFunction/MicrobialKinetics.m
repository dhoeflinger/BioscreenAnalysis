function [ lag_time, max_spec_growth_rate, max_od, delta_OD_max, median_od_max, doubling_time, note, goodness ] = MicrobialKinetics(OD_values, time_interval, incubation_time, threshold, model, double_hump)
%MicrobialKinetics -  For a specific dataset, attempts to determine a group
%of statistics on the growth which occurred.   Plots the dataset along the
%way, with the best fit regression.



%  time interval assumed to be half hour time blocks

%set to zero initially
max_spec_growth_rate = 0;
note = '';
lag_time = 0;
doubling_time = 0;

[max_od, max_location] = max(OD_values);

[lag_time, max_spec_growth_rate, median_od_max, min_od, goodness] = FindRegressionCurve(OD_values,time_interval, incubation_time, model, double_hump);


delta_OD_max = median_od_max - min(OD_values(1:4)); %max minus initial (min of first 4 elements to remove ouliers)

%use the max specific growth rate to calc the doubling time t_d = ln(2)/u
doubling_time = log(2) / max_spec_growth_rate;


absolute_no_growth_threshold = threshold;

relative_no_growth_threshold = threshold + min_od;

if (goodness.rsquare < 0.98)
   note = 'Bad r^2, Check Regression Plot'; 
end


if (max_od < absolute_no_growth_threshold && max_od < relative_no_growth_threshold)
    note = 'No Growth Detected, Check Plot';    
    lag_time = 0;
    max_spec_growth_rate = 0;
    doubling_time = 0;
    clf;
    return;
end

if (max_od < absolute_no_growth_threshold)
    note = 'Growth did not reach absolute threshold';
    lag_time = 0;
    max_spec_growth_rate = 0;
    doubling_time = 0;
    clf;
    return;
end

if (max_od < relative_no_growth_threshold)
    note = 'Growth did not reach relative threshold';
    lag_time = 0;
    max_spec_growth_rate = 0;
    doubling_time = 0;
    clf;
    return
end

if (strcmp(double_hump, 'double_hump') && lag_time == -2)
   note = 'No Double Hump Detected'; 
   return;
end
%report both max OD and median filtered max OD to excel 

lag_time_str = sprintf('lag time = %f, max growth = %f, doubling time = %f', lag_time, max_spec_growth_rate, doubling_time);

legend(lag_time_str, 'location', 'SouthOutside');
