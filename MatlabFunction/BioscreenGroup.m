function BioscreenGroup( directory )
%BIOSCREENGROUP Runs the Bioscreen Function on every excel file in the
%directory specified with default parameters
%   Default parameters are:
%          detecting the max timepoint automatically
%          growth threshold = 0.3
%          model = modlogistic

runnable_files = dir([directory '/' '*.xlsx']);

for i = 1:length(runnable_files)
    Bioscreen([directory '/' runnable_files(i).name]);
end

end

