function [lag_time, msgr, max_od, min_od, goodness] = FindRegressionCurve(OD_values, time_interval, incubation_time, model, double_hump)

if (nargin < 3) 
    model = 'modlogistic';
end;

timepoints = (0:size(OD_values)-1) * time_interval + incubation_time;




timepoints_med = timepoints;





OD_values_med = OD_values;
filter_width = 1;
for i=1+filter_width:size(OD_values)-filter_width
    sel = i - filter_width : i + filter_width;
    OD_values_med(i) = median(OD_values(sel)); 
    
end;

min_od = min(OD_values_med);
[max_od, max_index] = max(OD_values_med);


if (strcmp(double_hump, 'double_hump'))
    smoothing_range = 5;
    peaks = '';
    while (smoothing_range > 1 && length(peaks) < 1)
        smoothed_od_values =  smooth(OD_values, smoothing_range);
        [peaks, locs] = findpeaks(smoothed_od_values(1:min(max_index-2, length(smoothed_od_values))), 'MINPEAKDISTANCE', 3);
        j = 1;
        if (length(peaks) > 0)
            peaks_thresholded = [];
            locs_thresholded = [];
            for i = 1:length(peaks)
                if (peaks(i) < max_od * .85)
                    peaks_thresholded(j) = peaks(i);
                    locs_thresholded(j) = locs(i);
                    j = j + 1;
                end
            end
            peaks = peaks_thresholded;
            locs = locs_thresholded;
        end        
        smoothing_range = smoothing_range - 1;
    end
    if  (length(peaks)<1)
        lag_time = -1;
        msgr = -1;
        max_od = -1;
        min_od = -1;
        goodness.sse = -1;
        goodness.rsquare = -1;
        goodness.dfe = -1;
        goodness.adjrsquare=-1; 
        goodness.rmse =-1;
        return;       
    end

    [~, location_min] = min( smoothed_od_values(locs(length(locs)):max_index));
    location_min = location_min + locs(length(locs));

    OD_values_dh(1) = OD_values(1);
    timepoints_dh(1) = timepoints(1);

    OD_values_dh(2: length(OD_values)- location_min + 2) = OD_values(location_min: length(OD_values));
    timepoints_dh(2: length(OD_values)- location_min + 2) = timepoints(location_min: length(OD_values)); 
   
    OD_values = OD_values_dh;
    timepoints = timepoints_dh;
    max_index = max_index - location_min + 2;
end
timepoints_orig = timepoints;
OD_values_adj(size(OD_values)) = 0;
%log_OD_values(size(OD_values)) = 0;

for i=1:length(OD_values)
    if (i < max_index)
        OD_values_adj(i) = OD_values(i); 
    else
        OD_values_adj(i) = max_od;      
    end
%    log_OD_values(i) = log(OD_values_adj(i) - min_od + 0.01);
end;


max_od_adj = max(OD_values_adj);







%Appl. Environ. Microbiol. June 1990 vol. 56 no. 6 1875-1881

%and...


%0 Journal Article
%D 2000
%@ 0178-515X
%J Bioprocess Engineering
%V 23
%N 6
%R 10.1007/s004490000209
%T Development of mathematical models (Logistic, Gompertz and Richards models) describing the growth pattern of Pseudomonas putida (NICM 2174)
%U http://dx.doi.org/10.1007/s004490000209
%I Springer-Verlag
%8 2000-12-01
%A Annadurai, G.
%A Rajesh Babu, S.
%A Srinivasamoorthy, V. R.
%P 607-612
%G English

if (strcmpi(model, 'gompertz'))
      
    %Gompertz curve
    func = fittype('A * exp( -exp(-C*(x-B)))+D');

    [reg_curve, goodness] = fit(timepoints', OD_values_adj', func, 'Lower', [0 0 0 -1], 'Upper', [100 100 2 3], 'StartPoint', [0.5 1.5 0.2 0.1]);

    inflection_point = reg_curve.B;
    msgr = reg_curve.C;

    offset = .01;

    %is it ok to offset these values to make the log work properly?
    lag_time = inflection_point * time_interval - (log((reg_curve(inflection_point) + offset) - log(OD_values_adj(1) + offset)) / msgr);


elseif (strcmpi(model, 'modgompertz'))
    %modified Gompertz curve
    func = fittype('A * exp(-exp(((B * exp(1))/ A) * (C - x) + 1)) + D');

    [reg_curve, goodness] = fit(timepoints', OD_values_adj', func, 'Lower', [0 0.000001 0 0], 'Upper', [100 100 2 3], 'StartPoint', [0.5 1.5 0.2 0.1]);

    lag_time = reg_curve.B;
    msgr = reg_curve.C;


elseif (strcmpi(model, 'logistic'))
    %Logistic curve
    func = fittype('A / (1+exp(-C*(x-B))) + D');
    
    [reg_curve, goodness] = fit(timepoints', OD_values_adj', func, 'Lower', [0 0 0 0], 'Upper', [100 100 2 3], 'StartPoint', [0.5 1.5 0.2 0.1]);

    inflection_point = reg_curve.B;
    msgr = reg_curve.C;

    offset = .01;

    %is it ok to offset these values to make the log work properly?
    lag_time = inflection_point * time_interval - (log((reg_curve(inflection_point) + offset) - log(OD_values_adj(1) + offset)) / msgr);


    
elseif (strcmpi(model, 'modlogistic'))

    %Modified Logistic curve
    func = fittype('A / (1 + exp(((4*C)/A) * (B - x) + 2)) + D');
    
    [reg_curve, goodness] = fit(timepoints', OD_values_adj', func, 'Lower', [0 0.000001 0 0], 'Upper', [100 100 3 3], 'StartPoint', [2 1.5 0.2 0.1]);
    size_tmp = size(OD_values_adj);
%     for id=1:size_tmp
%         timepoints_t = timepoints;
%         OD_values_adj_t = OD_values_adj;
%         timepoints_t(id) = [];
%         OD_values_adj_t(id) = [];
%         [reg_curve_t, goodness_t] = fit(timepoints_t', OD_values_adj_t', func, 'Lower', [0 -0.5 0 0], 'Upper', [100 100 2 3], 'StartPoint', [0.5 1.5 0.2 0.1]);
%         if (goodness.rmse - goodness_t.rmse > 0.001)
%             timepoints = timepoints_t;
%             OD_values_adj = OD_values_adj_t;
%             goodness = goodness_t;
%             reg_curve = reg_curve_t;
%             size_tmp = size_tmp -1;
%             id
%             id = id -1;            
%         end
%             
%     end;
    
    lag_time = reg_curve.B;
    msgr = reg_curve.C;

    
%elseif (strcmpi(model, 'weibull'))
%    %weibull
%    func = fittype('B * exp(C * (log(x) - log(A)))');
%    
%    reg_curve = fit(timepoints', OD_values_adj', func, 'Lower', [0 0 0], 'Upper', [100 100 10], 'StartPoint', [0.5 1.5 0.5]);
%    
%    inflection_point = reg_curve.B;
%    msgr = reg_curve.C;
%
%    offset = .01;
%
%    %is it ok to offset these values to make the log work properly?
%    lag_time = inflection_point * time_interval - (log((reg_curve(inflection_point) + offset) - log(OD_values_adj(1) + offset)) / msgr);
%
end;

lag_time = max([lag_time 0]);

xlist = [lag_time lag_time];
ylist = [max_od_adj/4 0];

plot (reg_curve, timepoints_orig, OD_values);
hold on 
plot (timepoints_med, OD_values_med);
hold on 
%plot (timepoints, log_OD_values);
line(xlist, ylist);


hold off
