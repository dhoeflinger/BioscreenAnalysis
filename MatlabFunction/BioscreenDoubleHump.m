function BioscreenDoubleHump( excel_file_name, max_timepoint, growth_threshold, model)
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


Bioscreen( excel_file_name, max_timepoint, growth_threshold, model, 'double_hump');




end

