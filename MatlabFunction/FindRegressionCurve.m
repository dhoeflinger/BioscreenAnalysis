function [lag_time msgr] = FindRegressionCurve(OD_values, time_interval)
% size_od = size(OD_values);
% 
% for i=1:3
%     OD_values(size_od + i) = OD_values(size_od);
% end

timepoints = (1:size(OD_values)) * time_interval;
min_od = min(OD_values);


for i=1:size(OD_values)
    OD_values_adj(i) = OD_values(i) - min_od;
end

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


%Gompertz curve
%func = fittype('A * exp( -exp(B - C * x))');

%Logistic curve
A = max_od_adj;
func = fittype('A / (1+exp(-C*(x-B)))');


%weibull
%func = fittype('B * exp(C * (log(x) - log(A)))');


%lower and upper bounds are set up for logistic curve, 
%would change for gompertz
reg_curve = fit(timepoints', OD_values_adj', func, 'Lower', [0 1 0], 'Upper', [10 100 5], 'StartPoint', [0.5 1.5 0.5]);

fx = differentiate (reg_curve, timepoints);

[max_val max_idx] = max(fx);

inflection_point = reg_curve.B;
msgr = reg_curve.C;

offset = .1;

%is it ok to offset these values to make the log work properly?
lag_time = max_idx * time_interval - (log((reg_curve(inflection_point) + offset) - log(OD_values_adj(1) + offset)) / msgr);

xlist = [lag_time lag_time];
ylist = [reg_curve(inflection_point) OD_values_adj(1)];



plot (reg_curve, timepoints, OD_values_adj);

line(xlist, ylist);


