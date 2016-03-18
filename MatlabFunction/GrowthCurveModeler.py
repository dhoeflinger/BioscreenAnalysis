import os;
import glob;
import xlrd;
import xlsxwriter;
import numpy as np;
import scipy.optimize;
import scipy.stats;
import pylab as pl;
import matplotlib.pyplot as mpl;
import matplotlib.legend as mpll;
import warnings
from scipy.optimize import OptimizeWarning


def gompertz(x, A, B, C, D):
    return A * np.exp( - np.exp( -C * (x-B))) + D;    


def gompertz_range(x, A, B, C, D):
    if (A < 0  or A > 100):
        return 1e10;
    if (B < 0 or B > 100):
        return 1e10;
    if (C < 0 or C > 2):
        return 1e10;
    if (D < -1 or D > 3):
        return 1e10;
    return gompertz(x,A,B,C,D);


def modgompertz(x, A, B, C, D):
    return A * np.exp(-np.exp(((B * np.exp(1))/A) * (C-x) + 1)) + D;
    
def modgompertz_range(x, A, B, C, D):
    if (A < 0  or A > 100):
        return 1e10;
    if (B < 0.000001 or B > 100):
        return 1e10;
    if (C < 0 or C > 2):
        return 1e10;
    if (D < 0 or D > 3):
        return 1e10;    
    return modgompertz(x,A,B,C,D);


def logistic(x, A, B, C, D):
    return A / (1 + np.exp(-C * (x - B))) + D;

def logistic_range(x,A,B,C,D):
    if (A < 0  or A > 100):
        return 1e10;
    if (B < 0 or B > 100):
        return 1e10;
    if (C < 0 or C > 2):
        return 1e10;
    if (D < 0 or D > 3):
        return 1e10;
    return logistic(x,A,B,C,D);

def modlogistic(x,A,B,C,D):
    return A / (1 + np.exp(((4 * C)/A) * (B-x) + 2)) + D;

def modlogistic_range(x, A, B, C, D):
    # enforcing ranges on the values
    if (A < 0  or A > 100):
        A = A*3;
    if (B < 0.000001 or B > 100):
        B = B*3;
    if (C < 0 or C > 3):
       C = C*3;
    if (D < 0 or D > 3):
        D = D* 3;

    return modlogistic(x,A,B,C,D);        
    


def plot_results(full_plot, full_filename, lag_time=[], timepoints=[], time = [], regression=[], OD_values=[], msgr=[], doubling_time=[], delta_OD_max=[], min_od=[], median_od_max=[]):
    fig, ax = mpl.subplots();
            
    mpl.xlabel('Time (Hours)');
    mpl.ylabel('OD (600nm)');

    if (full_plot):
        lagtimestart = [lag_time, -10];
        lagtimestop = [lag_time, 10];
        #plot timepoints
        ax.plot(timepoints, OD_values[0:len(OD_values)-1], 'r.');
        ax.plot(time, regression, 'b-');      
        ax.plot(*zip(lagtimestart, lagtimestop), color='green');
        
        mpl.ylim((min(-0.2,min_od), max(1.6,median_od_max)));
        
    
        lag_time_str = "lag_time = %f\nmax growth = %f\ndoubling time = %f\nDelta OD Max = %f" % (lag_time, msgr, doubling_time, delta_OD_max);
    
          
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.95, 0.05, lag_time_str, transform=ax.transAxes, fontsize = 14, verticalalignment='bottom', horizontalalignment='right', bbox=props);    
    
    mpl.savefig(full_filename + ".png", dpi =600, format="png");
   
    mpl.close(fig);
    #mpl.show();    
    return;
    
def FindRegressionCurve(OD_values, time_interval, incubation_time, model, double_hump, threshold, full_filename):
    """
    FindRegressionCurve - using the parameters specified, attempts to fit a 
    regression best fit line of the specified model to the data supplied.
    """
    
    msgr = 0;
    note = ""
    lag_time = 0;
    goodness = 1;
    OD_values_med = np.empty([len(OD_values)]);
    timepoints = [];
    for i in range(len(OD_values)-1):
        timepoints.append(i * time_interval + incubation_time);

#    for i in range(len(OD_values)):
#        OD_values_med[i] = max(0, OD_values[i]);
        

    timepoints_med = timepoints;
     
    filter_width = 3;
    
    OD_values_med = medfilt(OD_values, filter_width);

    min_od = np.min(OD_values_med);
    
    median_od_max = np.max(OD_values_med);

    max_index = np.argmax(OD_values_med);

    double_hump_found = False;
    
    if (double_hump and max_index > 3):

        peaks_indices = scipy.signal.find_peaks_cwt(OD_values, np.arange(1,10));
        peaks_thresholded = [];
        locs_thresholded = [];
        
        if len(peaks_indices)>0:
            for i in range(len(peaks_indices)):
                if (OD_values[peaks_indices[i]] < median_od_max * 0.85):
                    peaks_thresholded.append(OD_values[peaks_indices[i]]);
                    locs_thresholded.append(peaks_indices[i]);

        if (length(peaks_thresholded) > 0):
            
            min_value = OD_values[peaks_thresholded[len(peaks_thresholded)-1]];
            for ele in OD_values[range(peaks_thresholded[len(peaks_thresholded)-1], len(OD_values) - 1 )]:
                min_value = min(ele, min_value);
            
            location_min = OD_values.index(min_value);
            
            OD_values_dh[0] = OD_values[0];
            timepoints_dh[0] = timepoints[0];
            
            for i in range(1, len(OD_values)-location_min +1): 
                OD_values_dh[i] = OD_values[i + location_min - 1];
                timepoints_dh[i] = timepoints[i + location_min -1];
            
            OD_values = OD_values_dh;
            timepoints = timepoints_dh;
            max_index = max_index - location_min + 2;
            double_hump_found = True;
    
    timepoints_orig = timepoints;
    
    OD_values_adj = [];
    for i in range(1, len(OD_values)-1):
       if (i < max_index):
           OD_values_adj.append(OD_values_med[i]); # should we use the median filtered data or just the raw data (changed to median 7/6/14) looks like it is still not median...
       else:
           OD_values_adj.append(median_od_max);        
    
    
    
    
    max_od_adj = max(OD_values_adj);
    
    
    """
    Appl. Environ. Microbiol. June 1990 vol. 56 no. 6 1875-1881

    and...

    0 Journal Article
    D 2000
    @ 0178-515X
    J Bioprocess Engineering
    V 23
    N 6
    R 10.1007/s004490000209
    T Development of mathematical models (Logistic, Gompertz and Richards models) describing the growth pattern of Pseudomonas putida (NICM 2174)
    U http://dx.doi.org/10.1007/s004490000209
    I Springer-Verlag
    8 2000-12-01
    A Annadurai, G.
    A Rajesh Babu, S.
    A Srinivasamoorthy, V. R.
    P 607-612
    G English
    """
    
    
    if (model == 'gompertz'):
        fitfunc_range = gompertz;
        fitfunc = gompertz;
        initial_guess =  [0.5, 1.5, 0.2, 0.1];
    elif (model == 'modgompertz'):
        fitfunc_range = modgompertz;
        fitfunc = modgompertz;
        initial_guess = [0.5, 1.5, 0.2, 0.1];
    elif (model == 'logistic'):
        fitfunc_range = logistic;
        fitfunc = logistic;
        initial_guess = [0.5, 1.5, 0.2, 0.1];
    elif (model == 'modlogistic'):
        fitfunc_range = modlogistic;
        fitfunc = modlogistic;
        initial_guess = [2, 1.5, 0.2, 0.1];
    else:
        print("Unsupported Model");
        return -1
        
    
    delta_OD_max = median_od_max - min(OD_values[0],OD_values[1], OD_values[2], OD_values[3]);
    
    if (delta_OD_max < threshold):
        note = 'No Growth Detected, Check Plot';    
        lag_time =0;
        msgr = 0;
        doubling_time = 0;
        noreg = 1;    
        goodness = 0;
        plot_results(False, full_filename);   
        return (lag_time, msgr, median_od_max, delta_OD_max, doubling_time, goodness, note);
        
    warnings.simplefilter("error", OptimizeWarning);
    try:
        (coef, est_cov) = scipy.optimize.curve_fit(fitfunc_range, timepoints[0:(len(timepoints)-1)], OD_values_adj, initial_guess, maxfev = 3000);
    except OptimizeWarning:
        print ('Maxed out cant estimate covariance, cannot fit curve');
        note = 'failed to fit curve';
        lag_time =0;
        msgr = 0;
        doubling_time = 0;
        noreg = 1;    
        goodness = 0;
        delta_OD_max = 0;
        plot_results(False, full_filename);   
        return (lag_time, msgr, median_od_max, delta_OD_max, doubling_time, goodness, note);
    except RuntimeError:        
        print ('Maxed out calls, cannot fit curve');
        note = 'failed to fit curve';
        lag_time =0;
        msgr = 0;
        doubling_time = 0;
        noreg = 1;    
        goodness = 0;
        delta_OD_max = 0;
        plot_results(False, full_filename);   
        return (lag_time, msgr, median_od_max, delta_OD_max, doubling_time, goodness, note);
    fit_values = fitfunc(timepoints[0:(len(timepoints)-1)],coef[0], coef[1],coef[2], coef[3]);    
    (slope, intercept, rsquared, pvalue, stderr) = scipy.stats.linregress(OD_values_adj, fit_values);
    if (model == 'gompertz' or model == 'logistic'):
        inflection_point = coef[1];
        msgr = coef[2];
        offset = 0.1;
        lag_time = inflection_point * time_interval - (np.log(fitfunc(inflection_point,coef[0],coef[1],coef[2], coef[3]) + offset) - np.log(OD_values_adj[0] + offset)) / msgr;              
    elif (model == 'modgompertz' or model == 'modlogistic'):
        lag_time = coef[1];
        msgr = coef[2];


        
    doubling_time = np.log(2) / msgr;
    
    lag_time = max(lag_time, 0);
    
    
 
    
    if (double_hump==1 and double_hump_found == 0):
        note = 'No Double Hump Detected'; 
        
    time = pl.arange(0, (len(OD_values)-1) * time_interval + incubation_time, 0.01);
    regression = modlogistic(time, coef[0],coef[1], coef[2], coef[3]);

    plot_results(True, full_filename, lag_time, timepoints, time, regression, OD_values, msgr, doubling_time, delta_OD_max, min_od, median_od_max);   
    
    #plot curve
    
    #plot line on lag_time
    
    #scale, and label
    
    return (lag_time, msgr, median_od_max, delta_OD_max, doubling_time, rsquared, note);
    

def medfilt (x, k):
    """Apply a length-k median filter to a 1D array x.
    Boundaries are extended by repeating endpoints.
    """
    assert k % 2 == 1, "Median filter length must be odd."
    assert x.ndim == 1, "Input must be one-dimensional."
    k2 = (k - 1) // 2
    y = np.zeros ((len (x), k), dtype=x.dtype)
    y[:,k2] = x
    for i in range (k2):
        j = k2 - i
        y[j:,i] = x[:-j]
        y[:j,i] = x[0]
        y[:-j,-(i+1)] = x[j:]
        y[-j:,-(i+1)] = x[-1]
    return np.median (y, axis=1)
    
    

def MicrobialKinetics(OD_values, time_interval, incubation_time, threshold, model, double_hump, full_filename):
    """
    MicrobialKinetics -  For a specific dataset, attempts to determine a group
    of statistics on the growth which occurred.   Plots the dataset along the 
    way, with the best fit regression.
    """
    # time interval assumed to be half hour time blocks


    #set to zero initially   
    max_od  = np.max(OD_values);    
    max_location = np.argmax(OD_values);
    
    (lag_time, max_spec_growth_rate,median_od_max, delta_OD_max, doubling_time, goodness, note) = FindRegressionCurve(
                    OD_values, time_interval, incubation_time, model, double_hump, threshold, full_filename);
    
    
    #report both max OD And median filtered max OF to excel


    
    
    #legend(lag_time_str, 'location', 'SouthOutside');
    
    
    return (lag_time, max_spec_growth_rate, max_od, median_od_max, delta_OD_max, doubling_time, goodness, note);
    
    
    
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
    print(file_or_dir);
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

    for k,v in varargin.items():
        if (k=='MaxTimepoint'):
            max_timepoint = v;
        if (k=='Threshold'):
            threshold = v;
        if (k=='Model'):
            model = v;
        if (k=='IncubationTime'):
            incubation_time = v;
        if (k=='DoubleHump'):
            double_hump = v;
    
    (path, file) = os.path.split(file_or_dir);
    (stub, ext) = os.path.splitext(file);
    
    if (path):
        path = path + "/";
        
    plots_folder = path + "results/" ;        
        
    if (not os.path.exists(plots_folder)):
        os.mkdir(plots_folder);
        
    plots_folder = plots_folder + stub + " plots/"

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

    output_file = path + 'results/' + stub + ' results.xlsx';
    output_workbook = xlsxwriter.Workbook(output_file);
    output_sheet = output_workbook.add_worksheet("Results");
    
    output_sheet.write_row(0,0, output);
        
    title_data0 = sheet.row(0);
    title_data1 = sheet.row(1);
    for i in range(num_columns-1):
        data_column = sheet.col_slice(i, 2);
        data_column_values = [];
        for d in data_column:
            data_column_values.append(d.value);
        if (sheet.cell(0,i)):
            if (not first_sugar):
                Strain_counts.append(strain_count);
            strain_count = 0;
            sugar_start_idx = i;
            first_sugar = False;
            sugar_count = sugar_count + 1;
            sugar = title_data0[i].value;
            Start_idxs.append(sugar_start_idx);
            Sugars.append(sugar)
        if (title_data1[i].value=="Time"):
            time_interval = data_column_values[1] - data_column_values[0];
        else:
            strain_count = strain_count + 1;
            strain = title_data1[i].value;
            sugar_folder = plots_folder + "/" + sugar;
            if (not os.path.exists(sugar_folder)):
                os.mkdir(sugar_folder);
            full_filename = sugar_folder + "/" +sugar + '-' + strain;
            
            if (max_timepoint < 0):
                (lagtime, max_u, OD_max, median_OD_max, delta_OD_max, doubling_time, goodness, note) = MicrobialKinetics(
                    np.array(data_column_values), time_interval, incubation_time, growth_threshold, model, double_hump, full_filename);
            else:
                (lagtime, max_u, OD_max, median_OD_max, delta_OD_max, doubling_time, goodness, note) = MicrobialKinetics(
                    np.array(data_column_values[0:max_timepoint / time_interval]), time_interval, incubation_time, growth_threshold, model, double_hump, full_filename);          
            lag_times.append(lagtime);
            output_sheet.write_row(i, 0, (sugar,strain,lagtime, max_u, doubling_time, OD_max,median_OD_max, delta_OD_max, \
                                          note, goodness));

            
            
            #plot title = name         
             
            #save plot to sugar folder
            #close plot
    
    
    
 

    


    