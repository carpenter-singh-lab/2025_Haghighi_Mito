import pandas as pd
from sqlalchemy import create_engine
from functools import reduce
import gc
import os
import time


def read_per_well_data(
    input_data_dir,
    annot,
    prof_workspace_folder_name="profiles",
    fformat=".parquet",
):
    batches = annot["Batch"].unique()

    df_agg_all_batches_ls = []
    for b in batches:
        print(b)
        #         if "Metadata_Source" in annot.columns:
        source_str = annot.loc[
            annot["Batch"] == b, "Metadata_Source"
        ].unique()[0]
        #             print(source_str)
        profile_path = (
            input_data_dir
            + source_str
            + "/workspace/"
            + prof_workspace_folder_name
            + "/"
        )
        #         else:
        #             profile_path = input_data_dir + "/workspace/profiles/"

        df_sag_ls = []
        plates_exist = os.listdir(profile_path + b)
        plates_meta = annot.loc[annot["Batch"] == b, "Metadata_Plate"].unique()
        plates = set(plates_meta) & set(plates_exist)
        for p in plates:
            print(p)

            fileName = profile_path + b + "/" + p + "/" + p + fformat
            #             print(fileName)
            if os.path.exists(fileName):
                if fformat == ".parquet":
                    sc_df = pd.read_parquet(fileName)
                elif fformat in [".csv", ".csv.gz"]:
                    sc_df = pd.read_csv(fileName)

                #         per_site_aggregate=sc_df.groupby(['Metadata_Well','Metadata_Site']).mean()[feature_list+['Count_Cells']].reset_index()
                sc_df["Metadata_Batch"] = b
                sc_df["Metadata_Plate"] = p
                df_sag_ls.append(sc_df)
                del sc_df
                gc.collect()
            else:
                print(fileName, " not exists")

        if df_sag_ls:
            df_sag = pd.concat(df_sag_ls, axis=0)
            df_agg_all_batches_ls.append(df_sag)

    df_agg_all_batches = pd.concat(
        df_agg_all_batches_ls, axis=0, ignore_index=True
    )
    return df_agg_all_batches


def read_per_well_data_csvs(input_data_dir,annot):
    batches=annot['Batch'].unique()
    
    df_agg_all_batches_ls=[]
    for b in batches:
        print(b)
        df_sag_ls=[]
        plates_exist=os.listdir(input_data_dir+b)
        plates_meta=annot.loc[annot['Batch']==b,'Metadata_Plate'].unique()
        plates=set(plates_meta) & set(plates_exist)
        for p in plates:
            print(p)
            
            fileName=input_data_dir+b+'/'+p+'/'+p+'.csv'
#             print(fileName)
            if os.path.exists(fileName):
                sc_df=pd.read_csv(fileName)

        #         per_site_aggregate=sc_df.groupby(['Metadata_Well','Metadata_Site']).mean()[feature_list+['Count_Cells']].reset_index()
                sc_df['Metadata_Batch']=b
                sc_df['Metadata_Plate']=p
                df_sag_ls.append(sc_df)
                del sc_df
                gc.collect()
            else:
                print(fileName,' not exists')

        if df_sag_ls:
            df_sag=pd.concat(df_sag_ls,axis=0)
            df_agg_all_batches_ls.append(df_sag)

    df_agg_all_batches=pd.concat(df_agg_all_batches_ls,axis=0)
    return df_agg_all_batches


def sample_single_cells_from_sql(input_data_dir,annot):
    batches=annot['Batch'].unique()
    
    df_agg_all_batches_ls=[]
    for b in batches:
        print(b)
        df_sag_ls=[]
        plates_exist=os.listdir(input_data_dir+b)
        plates_meta=annot.loc[annot['Batch']==b,'Metadata_Plate'].unique()
        plates=list(set(plates_meta) & set(plates_exist))
        for p in plates[:20]:
            
            fileName=input_data_dir+b+'/'+p+'/'+p+'.sqlite'
            print(p,fileName)
            n_rand_ims=100
            sc_df=read_single_cell_sql.readSingleCellData_sqlalch_random_image_subset(fileName,n_rand_ims)
    #         per_site_aggregate=sc_df.groupby(['Metadata_Well','Metadata_Site']).mean()[feature_list+['Count_Cells']].reset_index()
            sc_df['Metadata_Batch']=b
            sc_df['Metadata_Plate']=p
            df_sag_ls.append(sc_df)
            del sc_df
            gc.collect()

        df_sag=pd.concat(df_sag_ls,axis=0)
        df_agg_all_batches_ls.append(df_sag)

    df_agg_all_batches=pd.concat(df_agg_all_batches_ls,axis=0)
    return df_agg_all_batches


def readSingleCellData_sqlalch_features_subset(fileName, feature_list):
    start1 = time.time()

    d = {}
    for f in feature_list:
        comp = f.split("_")[0].lower()
        if comp in d.keys():
            d[comp] += [f]
        else:
            d[comp] = [f]

    sql_file = "sqlite:////" + fileName
    engine = create_engine(sql_file)
    conn = engine.connect()

    compartments = list(d.keys())

    plateDf_list = []
    for compartment in compartments:
        features = ", ".join(d[compartment])
        query_cols = "TableNumber, ImageNumber, ObjectNumber, " + features  # +", "+f2
        compartment_query = "select {} from {}".format(query_cols, compartment)
        plateDf_list.append(pd.read_sql(sql=compartment_query, con=conn))

    plateDf = reduce(
        lambda left, right: pd.merge(
            left, right, on=["TableNumber", "ImageNumber", "ObjectNumber"]
        ),
        plateDf_list,
    )
    #     plateDf=plateDf.dropna()
    img_query = "select * from {}".format("Image")
    plateImageDf = pd.read_sql(sql=img_query, con=conn)
    #     print(plateImageDf.columns)
    plateDfwMeta = pd.merge(plateDf, plateImageDf, on=["TableNumber", "ImageNumber"])

    end1 = time.time()
    print("time elapsed:", (end1 - start1) / 60, " mins")
    return plateDfwMeta


def form_per_site_aggregated_profiles(annot,input_data_dir,output_dir,feature_list2):
    import pandas as pd
    from sqlalchemy import create_engine
    from functools import reduce
    import gc

    batches=annot['Batch'].unique().tolist()
    for b in batches:
        
        df_sag_ls=[]
        if "Metadata_Source" in annot.columns:
            src=annot.loc[annot['Batch']==b,'Metadata_Source'].unique()[0]
            print(b,src)
            input_data_dir = '/'.join([src if 'source_' in i else i for i in input_data_dir.split('/')])
            
        plates_exist=os.listdir(input_data_dir+b)
        plates_meta=annot.loc[annot['Batch']==b,'Metadata_Plate'].unique()
        plates=set(plates_meta) & set(plates_exist)
        for p in plates:
            print(p)
            fileName=input_data_dir+b+'/'+p+'/'+p+'.sqlite'
            print(fileName)
            sc_df=readSingleCellData_sqlalch_features_subset(fileName,feature_list2)
            cell_count_col_name=sc_df.columns[sc_df.columns.str.contains('Count_Cell')].values[0]
            per_site_aggregate=sc_df.groupby(['Metadata_Well','Metadata_Site']).mean(numeric_only=True)[feature_list2+\
                                                                            [cell_count_col_name]].reset_index()
#                                                                         ['Count_Cells']].reset_index()
#                                                                   ['Count_CellsIncludingEdges']].reset_index() for crispr
                                                                                       
            per_site_aggregate['Count_Cells']=per_site_aggregate[cell_count_col_name]                                                                      
            per_site_aggregate['Metadata_Batch']=b
            per_site_aggregate['Metadata_Plate']=p
            df_sag_ls.append(per_site_aggregate)
#             del sc_df
#             gc.collect()

        df_sag=pd.concat(df_sag_ls,axis=0)
        fileNameToSave=output_dir+'/'+b+"_site_agg_profiles"
        print(fileNameToSave)
        saveDF_to_CSV_GZ_no_timestamp(df_sag,fileNameToSave)
        
    return


def saveDF_to_CSV_GZ_no_timestamp(df,filename):
    from gzip import GzipFile
    from io import TextIOWrapper
    with TextIOWrapper(GzipFile(filename+'.csv.gz', 'w', mtime=0), encoding='utf-8') as fd:
        df.to_csv(fd,index=False,compression='gzip')
        
    return


def saveAsNewSheetToExistingFile(filename, newDFs, newSheetNames, keep_index_column=True):
    # Ensure newDFs and newSheetNames are lists
    if not isinstance(newDFs, list):
        newDFs = [newDFs]
    if not isinstance(newSheetNames, list):
        newSheetNames = [newSheetNames]
    
    # Check that the number of DataFrames matches the number of sheet names
    if len(newDFs) != len(newSheetNames):
        raise ValueError("The number of DataFrames must match the number of sheet names.")
    
    if os.path.exists(filename):
        excel_book = pxl.load_workbook(filename)

        with pd.ExcelWriter(filename, engine='openpyxl', mode='a') as writer:
            writer.book = excel_book

            # Write each DataFrame to its corresponding sheet
            for df, sheet_name in zip(newDFs, newSheetNames):
                df.to_excel(writer, sheet_name=sheet_name, index=keep_index_column)
    else:
        with pd.ExcelWriter(filename, engine='openpyxl') as writer:
            for df, sheet_name in zip(newDFs, newSheetNames):
                df.to_excel(writer, sheet_name=sheet_name, index=keep_index_column)
        
    return

from scipy.stats import f
def TwoSampleT2Test(X, Y):
    nx, p = X.shape
    ny, _ = Y.shape
    delta = np.mean(X, axis=0) - np.mean(Y, axis=0)
    Sx = np.cov(X, rowvar=False)
    Sy = np.cov(Y, rowvar=False)
    S_pooled = ((nx-1)*Sx + (ny-1)*Sy)/(nx+ny-2)
    S_pooled = S_pooled + np.eye(S_pooled.shape[0]) * 1e-6
    t_squared = (nx*ny)/(nx+ny) * np.matmul(np.matmul(delta.transpose(), np.linalg.inv(S_pooled)), delta)
    statistic = t_squared * (nx+ny-p-1)/(p*(nx+ny-2))
    F = f(p, nx+ny-p-1)
    p_value = 1 - F.cdf(statistic)
#     print(f"Test statistic: {statistic}\nDegrees of freedom: {p} and {nx+ny-p-1}\np-value: {p_value}")

    # Convert F-statistic to z-score
    z_score = (statistic - (p / (nx + ny - p - 1))) / np.sqrt((2 * p * (nx + ny - p - 1)) / ((nx + ny - 2) * (nx + ny - p - 1)))
    std_p_val = 2 * (1 - norm.cdf(abs(z_score)))
    
    return statistic, p_value, std_p_val

import numpy as np
from scipy.stats import f, chi2

def HotellingsT_internal(X, Y, test='f'):
    n1, p = X.shape
    n2 = Y.shape[0]
    
    mu=np.zeros(p)
    
    # Calculate means and differences
    Xmeans = np.mean(X, axis=0)
    Ymeans = np.mean(Y, axis=0)
    X_diff = X - Xmeans
    Y_diff = Y - Ymeans
    
    # Calculate pooled covariance matrix
    S_pooled = 1 / (n1 + n2 - 2) * (X_diff.T @ X_diff + Y_diff.T @ Y_diff)
    
    # Calculate test statistic
    diff_means = Xmeans - Ymeans - mu
    if test == 'f':
        test_statistic = n1 * n2 / (n1 + n2) * diff_means @ np.linalg.inv(S_pooled) @ diff_means.T * (n1 + n2 - p - 1) / (p * (n1 + n2 - 2))
        df1 = p
        df2 = n1 + n2 - p - 1
        p_value = 1 - f.cdf(test_statistic, df1, df2)
    elif test == 'chi':
        test_statistic = n1 * n2 / (n1 + n2) * diff_means @ np.linalg.inv(S_pooled) @ diff_means.T
        df1 = p
        df2 = None
        p_value = 1 - chi2.cdf(test_statistic, df1)
    else:
        return "Invalid test type"
    
    return test_statistic,  p_value

    
    
from scipy.signal import find_peaks
def find_end_slope(data, height=None):
    peaks, _ = find_peaks(data, height=height,width=2)
    valleys, _ = find_peaks(-data, height=height,width=2)
    extermas=np.concatenate((peaks, valleys))
    if extermas.size==0:
        return np.nan,np.nan
    
    last_peak_ind=np.max(extermas)
    slope=data[-1]-data[last_peak_ind]
    return last_peak_ind,slope




import matplotlib.pyplot as plt

from scipy.signal import find_peaks

from scipy.signal import savgol_filter


def smooth_data(data, window_length=5, polyorder=3):
    return savgol_filter(data, window_length, polyorder)


def find_end_slope(data, height=None):
    peaks, _ = find_peaks(data, height=height, width=1)
    valleys, _ = find_peaks(-data, height=height, width=1)
    extermas = np.concatenate((peaks, valleys))
    if extermas.size == 0:
        return np.nan, 0

    last_peak_ind = np.max(extermas)
    slope = data[-1] - data[last_peak_ind]
    return last_peak_ind, slope


def subtract_control(group, control_df_perplate):
    batch_plate = group.name
    control_values = control_df_perplate.loc[batch_plate]
    return group - control_values


def find_end_slope2(data, height=None):
    data = smooth_data(data)
    #     min_max_indc = [np.argmax(data[3:] + 3), np.argmin(data[3:] + 3)]
    min_max_indc = [np.argmax(data), np.argmin(data)]
    last_peak_ind0 = [i for i in min_max_indc if i < len(data) - 2]
    if last_peak_ind0 == []:
        return 0, 0
    last_peak_ind = np.max(last_peak_ind0)
    last_two_points_amplitude = (data[-1] + data[-2]) / 2
    slope = (last_two_points_amplitude - data[last_peak_ind]) / (
        len(data) - last_peak_ind - 1
    )
    return last_peak_ind, slope





################################################################################
def handle_nans(
    df_input, cp_features, thrsh_null_ratio=0.05, thrsh_std=0.0001, fill_na_method=None
):
    """
    from the all df_input columns extract cell painting measurments
    the measurments that should be used for analysis

    Inputs:
    df_input: dataframes with all the annotations available in the raw data
    fill_na_method: 'interpolate' or 'median' or None
                    interpolate makes sense for single cell data of arrayed experiment since it fills NA values
                    to the nearest cell values
    Outputs: cp_features, cp_features_analysis

    """

    import pandas as pd
    import numpy as np

    #     cp_features=df_input.columns[df_input.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")].tolist()

    print("cp_features:", len(cp_features))
    object_type_columns = df_input[cp_features].select_dtypes([object]).columns
    df_input[object_type_columns] = df_input[object_type_columns].apply(
        pd.to_numeric, errors="coerce"
    )

    #     df_input=df_input.replace([np.inf, -np.inf,'nan'], np.nan)
    #     df_input[cp_features] = df_input[cp_features].apply(pd.to_numeric, errors='coerce')
    #     df_input[cp_features] = df_input[cp_features].astype(float)

    #     thrsh_null_ratio=0.05; thrsh_std=0.0001;
    cols2remove_manyNulls = [
        i
        for i in cp_features
        if (df_input[i].isnull().sum(axis=0) / df_input.shape[0]) > thrsh_null_ratio
    ]
    cols2remove_lowVars = (
        df_input[cp_features]
        .std()[df_input[cp_features].std() < thrsh_std]
        .index.tolist()
    )

    cols2removeCP = cols2remove_manyNulls + cols2remove_lowVars
    print("cols2remove_manyNulls", cols2remove_manyNulls)

    print("cols2remove_lowVars", cols2remove_lowVars)

    cp_features_analysis = list(set(cp_features) - set(cols2removeCP))
    print(
        "len cp_features_analysis/nan cols/low vars:",
        len(cp_features_analysis),
        len(cols2remove_manyNulls),
        len(cols2remove_lowVars),
    )
    #     cp_features_analysis_filt1 = list(set(cp_features) - set(cols2removeCP))
    #     df_numeric_columns = df_input.select_dtypes([np.number]).columns
    #     cp_features_analysis = list(set(df_numeric_columns) & set(cp_features_analysis_filt1))

    df_p_s = df_input.drop(cols2removeCP, axis=1)

    #     print(cols2removeCP)

    #     df_p_s[cp_features_analysis] = df_p_s[cp_features_analysis].interpolate()
    if fill_na_method == "median":
        df_p_s.loc[:, cp_features_analysis] = df_p_s.loc[
            :, cp_features_analysis
        ].fillna(df_p_s[cp_features_analysis].median())
    elif fill_na_method == "interpolate":
        df_p_s.loc[:, cp_features_analysis] = df_p_s.loc[
            :, cp_features_analysis
        ].interpolate()
    elif fill_na_method == "drop-rows":
        print("before dropping nan rows: ", df_p_s.shape)
        #         print('1',df_p_s.shape,df_p_s.dropna(subset=cp_features_analysis).reset_index(drop=True).shape)
        df_p_s = df_p_s.dropna(subset=cp_features_analysis).reset_index(drop=True)
        print("after dropping nan rows: ", df_p_s.shape)

    elif (
        fill_na_method == "interpolate_sim_col"
    ):  # interpolate based on the columns with highest correlation
        print("Not implemented yet! Nothing got dropped! ")

    #     row_has_NaN = df_p_s[cp_features_analysis].isnull().any(axis=1)
    #     print(row_has_NaN)
    #     print(df_p_s[cp_features_analysis].dropna().shape,df_p_s[cp_features_analysis].shape)
    #     df_p_s[cp_features_analysis] = df_p_s[cp_features_analysis].dropna()
    #     dataframe.fillna(0)

    return df_p_s, cp_features_analysis


def standardize_per_catX(df, column_name, cp_features):
    # """

    df_scaled_perPlate = df.copy()
    group_means = df.groupby(column_name)[cp_features].mean()
    group_stds = df.groupby(column_name)[cp_features].std()
    df_scaled_perPlate[cp_features] = (
        df[cp_features] - group_means.loc[df[column_name]].values
    ) / group_stds.loc[df[column_name]].values
    return df_scaled_perPlate


##############################################
from scipy.stats import ttest_ind, norm
def cohens_d(x, y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    pooled_std = np.sqrt(((nx - 1) * np.std(x, ddof=1) ** 2 + (ny - 1) * np.std(y, ddof=1) ** 2) / dof)
    return (np.mean(x) - np.mean(y)) / pooled_std

# Function to convert t-statistic to z-score
def t_to_z(t_stat, df):
    return t_stat / np.sqrt(df / (df + t_stat**2))

# Function to calculate standardized p-value from z-score
def z_to_p(z):
    return 2 * (1 - norm.cdf(abs(z)))


def bh_adjusted_critical_value(pvalues, fdr=0.05):
    sorted_pvalues = np.sort(pvalues)
    m = len(pvalues)
    ranks = np.arange(1, m + 1)
    critical_values = (ranks / m) * fdr
    below_threshold = sorted_pvalues <= critical_values
    if np.any(below_threshold):
        adjusted_critical = sorted_pvalues[below_threshold].max()
    else:
        adjusted_critical = np.nan  # nan if no p-values are below the threshold
    return adjusted_critical