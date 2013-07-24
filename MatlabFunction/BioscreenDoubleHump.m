function BioscreenDoubleHump( excel_file_name, max_timepoint, growth_threshold, model, incubation_time)
%BIOSCREENDOUBLEHUMP Summary of this function goes here
%   Detailed explanation goes here


if (nargin < 2)
    max_timepoint = -1;
end

if (nargin < 3)
    growth_threshold = 0.3;
end

if (nargin < 4)
    model = 'modlogistic';
end

if (nargin < 5)
    incubation_time = 1.0;
end


Bioscreen( excel_file_name, max_timepoint, growth_threshold, model, incubation_time, 1);




end

