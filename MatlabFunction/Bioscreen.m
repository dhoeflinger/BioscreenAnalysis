function Bioscreen( excel_file_name, max_timepoint, growth_threshold, model)
%Bioscreen - Calculates some metrics for growth curves, as well as graphing
% a regression curve of the data points. 
%PARAMS: 
% excel_file_name - input data file, must be formatted properly
% max_timepoint - final timepoint of the dataset
% growth_threshold - threshold which describes the minimum OD reading to signify growth (default 0.2)
% model - regression model to plot the curves, options are 'modgompertz' (default),
%         'gompertz', 'logistic', and 'modlogistic'
% 

if (nargin < 2)
    max_timepoint = -1;
end

if (nargin < 3)
    growth_threshold = 0.2;
end

if (nargin < 4)
    model = 'modgompertz';
end

[Data, title_data] = xlsread(excel_file_name);
dims = size(title_data);

plots_folder = [ excel_file_name(1:(strfind(excel_file_name, '.')-1)) ' plots/'];
        
if ~exist(plots_folder, 'dir')
    mkdir(plots_folder);
end

output(1,:) = {'Sugar', 'Strain', 'Lag Time (hours)', 'Max Specific Growth Rate (1/hours)',  'Doubling Time (hours)', 'Max OD', 'Median OD', 'Notes', 'SSE' , 'R^2', 'DFE', 'ADJ R^2', 'RMSE'};
time_interval =0.5;

sugar = '';
strain_count = 0;
sugar_count = 0;

first_sugar = 1;
Sugars = {};
Start_idxs = [];
Strain_counts = [];
lag_times = [];
for i=1:dims(2)
    if (~isempty(char(title_data(1,i))))
        if  (~first_sugar)
            Strain_counts(sugar_count) = strain_count;
        end
        strain_count = 0;
        sugar_start_idx = i;
        first_sugar = 0;
        sugar_count = sugar_count + 1;
        
        sugar = char(title_data(1,i));
        Start_idxs(sugar_count) = sugar_start_idx;
        Sugars(sugar_count) = {char(sugar)};

    end
    if (strcmpi(char(title_data(2,i)),'Time'))
        time_interval = Data(2,i) - Data(1,i);
    else
        strain_count = strain_count + 1;
        strain = char(title_data(2,i));
        h = figure;
        if (max_timepoint < 0)
            [lagtime, max_u, OD_max, median_OD_max, doubling_time, note, goodness] = MicrobialKinetics(Data(:,i), time_interval, growth_threshold, model);
        else
            [lagtime, max_u, OD_max, median_OD_max, doubling_time, note, goodness] = MicrobialKinetics(Data(1:max_timepoint/time_interval,i), time_interval, growth_threshold, model);
        end
        lag_times(i) = lagtime;
        output(i,:) = {sugar,strain, lagtime, max_u, doubling_time, OD_max, median_OD_max, note, goodness.sse, goodness.rsquare, goodness.dfe, goodness.adjrsquare, goodness.rmse};
        name = [sugar '-' strain];
        title( name );
        sugar_folder = [plots_folder '/' sugar];
        if ~exist(sugar_folder, 'dir')
            mkdir(sugar_folder);
        end
        saveas(h,[sugar_folder '/' name], 'bmp');
        close(h);
    end
    
end

Strain_counts(sugar_count) = strain_count;


% for j = 1:sugar_count
%     name = char(Sugars(j));
%     start_idx = Start_idxs(j);
%     count = Strain_counts(j);
%     for k = start_idx:start_idx + count-1
%         figure;
%         if (max_timepoint < 0)
%             plot (Data(:,start_idx -1), Data(:,k),lag_times(k),0);
%         else
%             plot (Data(1:max_timepoint/time_interval,start_idx -1), Data(1:max_timepoint / time_interval,k), lag_times(k), 0);
%         end
%         xlabel('hours');
%         ylim([-0.5 2]);
%         title([name ' ' char(title_data(2,k))]);
%     end
% end

output_file = ['result ' excel_file_name]; 
xlswrite(output_file, output);



end

