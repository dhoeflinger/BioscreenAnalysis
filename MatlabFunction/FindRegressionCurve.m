function [lag_time, msgr, max_od, min_od, goodness] = FindRegressionCurve(OD_values, time_interval, model)

if (nargin < 3) 
    model = 'modgompertz';
end;

timepoints = (1:size(OD_values)) * time_interval;

timepoints_orig = timepoints;


OD_values_adj(size(OD_values)) = 0;
%log_OD_values(size(OD_values)) = 0;

OD_values_med = OD_values;
filter_width = 1;
for i=1+filter_width:size(OD_values)-filter_width
    sel = i - filter_width : i + filter_width;
    OD_values_med(i) = median(OD_values(sel)); 
    
end;

min_od = min(OD_values_med);
[max_od, max_index] = max(OD_values_med);

for i=1:size(OD_values_med)
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

    [reg_curve, goodness] = fit(timepoints', OD_values_adj', func, 'Lower', [0 -0.5 0 0], 'Upper', [100 100 2 3], 'StartPoint', [0.5 1.5 0.2 0.1]);

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
    
    [reg_curve, goodness] = fit(timepoints', OD_values_adj', func, 'Lower', [0 -0.5 0 0], 'Upper', [100 100 3 3], 'StartPoint', [0.5 1.5 0.2 0.1]);
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
plot (timepoints_orig, OD_values_med);
hold on 
%plot (timepoints, log_OD_values);
line(xlist, ylist);


hold off
