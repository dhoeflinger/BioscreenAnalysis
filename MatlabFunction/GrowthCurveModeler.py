import os;
import glob;
import xlrd:
import xlsxwriter;

def GrowthCurveModeler( file_or_dir, **varargin):
    """
    GrowthCurveModeler - Calculates some metrics for growth curves, as well as graphing
        a regression curve of the data points. 

    Required parameter: 
    
    file_or_dir - input data file or directory of input files, must be formatted properly

    Optional parameters:

    MaxTimepoint - value {default is total specified in dataset}
      Final timepoint of the dataset 

    Threshold - value {default 0.3}
      Threshold which describes the minimum OD reading to signify growth (default 0.3)
 
     Model - [{'modlogistic'} | 'gompertz' | 'logistic' | 'modgompertz']
       Regression model to plot the curves

     IncubationTime - value {default 1.0}
       If you incubate your cells before recording timepoints,
       this will appropriately shift your data points in time
       for lag time calculations 

     DoubleHump - [True | {False}] 
       Flag to indicate double hump processing, this expects all datasets to include a
       double/multi hump and should remove all growth humps before the "main curve" 
 


      Examples:
          GrowthCurveModeler('dataset.xlsx', DoubleHump=True, Threshold=0.2);

         GrowthCurveModeler('.', IncubationTime=1.5);
 
         GrowthCurveModeler('folder_containing_xlsxfiles');
    """

    if (os.path.isdir(file_or_dir)):
        runnable_files = glob.glob(file_or_dir + "/*.xlsx");
        for r in runnable_files:
            GrowthCurveModeler(r, varargin);
        return;

    max_timepoint = -1;
    growth_threshold = 0.3;
    model = 'modlogistic';
    incubation_time = 1.0;
    double_hump = 0; 

    for k,v in varargin:
        if (k="MaxTimepoint"):
            max_timepoint = v;
        if (k="Threshold"):
            threshold = v;
        if (k="Model"):
            model = v;
        if (k="IncubationTime"):
            incubation_time = v;
        if (k="DoubleHump"):
            double_hump = v;
    
    (path, file) = os.path.split(file_or_dir);
    (stub, ext) = os.path.splitext(file);
    
    if (path):
        path = path + "/";
        
    plots_folder = path + "results/" + stub + " plots/";        
        
    if (not os.path.exists(plots_folder)):
        os.mkdir(plots_folder);
    
    workbook = xlrd.open_workbook(file_or_dir);
    sheet = workbook.sheet_by_index(0);
    num_columns = sheet.ncols;    
    
    output = ('Sugar', 'Strain', 'Lag Time (hours)', 'Max Specific Growth Rate (1/hours)',\
            'Doubling Time (hours)', 'Max OD', 'Median OD', 'Delta OD', 'Notes', 'R^2', 'SSE' , 'RMSE');   
    
    time_interval = 0.5;
    sugar = '';
    strain_count = 0;
    sugar_count = 0;    
    
    first_sugar = True;
    Sugars = []
    Start_idxs = [];
    Strain_counts = [];
    lag_times = [];

    output_workbook = xlsxwriter.Workbook(output_file);
    output_sheet = output_workbook.add_worksheet("Results");
    
    output_sheet.write_row(0,0, output);
        
    title_data0 = sheet.row(0);
    title_data1 = sheet.row(1);
    for i = 0:num_columns-1:
        data_column = sheet.col_slice(i, 2);
        if (sheet.cell(0,i)):
            if (first_sugar):
                Strain_counts[sugar_count] 
            strain_count = 0;
            sugar_Start_idx = i;
            first_sugar = False;
            sugar_count = sugar_count + 1;
            sugar = title_data0[i];
            Start_idxs[sugar_count] = sugar_start_idx;
            Sugars[sugar_count] = sugar
        if (title_data1[i]=="Time"):
            time_interval = data_column[1] - data_column[0];
        else:
            strain_count = strain_count + 1;
            strain = title_data1(0);
            if (max_timepoint < 0):
                (lagtime, max_u, OD_max, median_OD_max, delta_OD_max, doubling_time, note, goodness) = MicrobialKinetics(
                    data_column, time_interval, incubation_time, growth_threshold, model, double_hump);
            else:
                (lagtime, max_u, OD_max, median_OD_max, delta_OD_max, doubling_time, note, goodness) = MicrobialKinetics(
                    data_column[0:max_timepoint / time_interval], time_interval, incubation_time, growth_threshold, model, double_hump);          
            lag_times[i] = lagtime;
            output_sheet.write_row(i, 0, (sugar,strain,lagtime, max_u, doubling_time, OD_max,median_OD_max, delta_OD_max, \
                                          note, goodnesss.rsquare, goodness.sse, goodness.rmse));
            name = sugar + '-' + strain;
            
            
            #plot title = name         
             
            sugar_folder = plots_folder + "/" + sugar;
            if (not os.path.exists(sugar)):
                os.mkdir(sugar_folder);
            
            #save plot to sugar folder
            #close plot
    
    
    
    
    
def MicrobialKenetics(OD_values, time_interval, incubation_time, threshold, model, double_hump):
    