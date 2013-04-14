function Bioscreen( excel_file_name, excel_output_file, max_timepoint, growth_threshold)
%READEXCEL Summary of this function goes here
%   Detailed explanation goes here

if (nargin < 3)
    max_timepoint = -1;
end

[Data, title_data] = xlsread(excel_file_name);
dims = size(title_data);

if ~exist('plots', 'dir')
    mkdir('plots');
end

output(1,:) = {'Sugar', 'Strain', 'Lag Time (hours)', 'Max Specific Growth Rate (1/hours)',  'Max OD', 'Doubling Time (hours)', 'Notes'};
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
            [lagtime, max_u, OD_max, doubling_time, note] = MicrobialKinetics(Data(:,i), time_interval, growth_threshold);
        else
            [lagtime, max_u, OD_max, doubling_time, note] = MicrobialKinetics(Data(1:max_timepoint/time_interval,i), time_interval, growth_threshold);
        end
        lag_times(i) = lagtime;
        output(i,:) = {sugar,strain, lagtime, max_u, OD_max, doubling_time, note};
        name = [sugar '-' strain];
        title( name );
        
        saveas(h,['plots/' name], 'bmp');
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

xlswrite(excel_output_file, output);



end

