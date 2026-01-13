import numpy as np
import pickle
import pandas as pd
from functools import reduce

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




################################################################################
def extract_cpfeature_names(df_input):
    """
    from the all df_input columns extract cell painting measurments
    the measurments that should be used for analysis

    Inputs:
    df_input: dataframes with all the annotations available in the raw data

    Outputs: cp_features, cp_features_analysis

    """

    cp_features = df_input.columns[
        df_input.columns.str.contains("Cells_|Cytoplasm_|Nuclei_")
    ].tolist()
    locFeature2beremoved = list(
        filter(
            lambda x: "_X" in x
            or "_Y" in x
            or "_Z" in x
            or "_x" in x
            or "_y" in x
            or "_z" in x,
            cp_features,
        )
    )
    metadataFeature2beremoved = list(filter(lambda x: "etadata" in x, cp_features))

    blackListFeatures = df_input.columns[
        df_input.columns.str.contains(
            "Nuclei_Correlation_Manders_"
            "|Nuclei_Correlation_RWC_|Nuclei_Granularity_14_|Nuclei_Granularity_15_|Nuclei_Granularity_16_"
        )
    ].tolist()

    if 0:  # changed to the above approach as the below fixed hard coded one was
        with open("./blackListFeatures.pkl", "rb") as f:
            blackListFeatures = pickle.load(f)

    cp_features_analysis = list(
        set(cp_features)
        - set(locFeature2beremoved)
        - set(metadataFeature2beremoved)
        - set(blackListFeatures)
    )

    return cp_features, cp_features_analysis


################################################################################
def find_correlation(data, threshold=0.9, remove_negative=False):

    """
    Inputs
    data : pandas DataFrame
    threshold : float
        correlation threshold, will remove one of pairs of features with a
        correlation greater than this value.
    remove_negative: Boolean
        If true then features which are highly negatively correlated will
        also be returned for removal.

    Output
    to_drop(list): list of column names to be removed
    """
    corr_mat = data.corr()
    if remove_negative:
        corr_mat = corr_mat.abs()

    upper = corr_mat.where(np.triu(np.ones(corr_mat.shape), k=1).astype(np.bool))

    to_drop = [column for column in upper.columns if any(upper[column] > threshold)]

    return to_drop


################################################################################
def readSingleCellData_sqlalch(fileName, compartments):
    from sqlalchemy import create_engine

    sql_file = "sqlite:////" + fileName
    engine = create_engine(sql_file)
    conn = engine.connect()
    #     compartments=["cells", "cytoplasm", "nuclei"]
    # compartments=["Neurites","CellBodies","CellBodiesPlusNeurites","Nuclei","Cytoplasm"]
    plateDf_list = []
    for compartment in compartments:
        compartment_query = "select * from {}".format(compartment)
        plateDf_list.append(pd.read_sql(sql=compartment_query, con=conn))

    plateDf = reduce(
        lambda left, right: pd.merge(
            left, right, on=["TableNumber", "ImageNumber", "ObjectNumber"]
        ),
        plateDf_list,
    )

    compartment_query = "select * from {}".format("Image")
    plateImageDf = pd.read_sql(sql=compartment_query, con=conn)

    plateDfwMeta = pd.merge(plateDf, plateImageDf, on=["TableNumber", "ImageNumber"])
    plateDfwMeta = plateDfwMeta.loc[:, ~plateDfwMeta.columns.duplicated()]

    return plateDfwMeta


################################################################################
def check_feature_similarity_dendrogram(data, feature_names, figsize):
    from scipy.cluster import hierarchy
    from scipy.spatial.distance import squareform
    import matplotlib.pyplot as plt

    #     import hdmedians as hd

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    #     corr = spearmanr(X).correlation
    corr = data[feature_names].corr().values

    # Ensure the correlation matrix is symmetric
    corr = (corr + corr.T) / 2
    np.fill_diagonal(corr, 1)

    # We convert the correlation matrix to a distance matrix before performing
    # hierarchical clustering using Ward's linkage.
    distance_matrix = 1 - np.abs(corr)
    dist_linkage = hierarchy.ward(squareform(distance_matrix))
    dendro = hierarchy.dendrogram(
        dist_linkage, labels=feature_names, ax=ax1, leaf_rotation=90
    )
    dendro_idx = np.arange(0, len(dendro["ivl"]))

    pos = ax2.imshow(corr[dendro["leaves"], :][:, dendro["leaves"]], vmin=-1, vmax=1)
    ax2.grid(False)
    fig.colorbar(pos, ax=ax2)
    ax2.set_xticks(dendro_idx)
    ax2.set_yticks(dendro_idx)
    ax2.set_xticklabels(dendro["ivl"], rotation="vertical")
    ax2.set_yticklabels(dendro["ivl"])
    fig.tight_layout()

    return fig


from typing import Optional
import matplotlib.pyplot as plt

def plot_lollipop_feature_corr(
    df: pd.DataFrame,
    target_feature: str,
    method: str = "spearman",
    top_n: Optional[int] = None,
    abs_corr: bool = False,
    ax: Optional[plt.Axes] = None,
    decimals: int = 2,
    label_fontsize: int = 8,
):
    """
    Lollipop chart of correlation between a target feature and all other features,
    with correlation values annotated next to each point.

    Text is placed on the right of positive correlations and on the
    left of negative correlations.
    """

    if target_feature not in df.columns:
        raise ValueError(f"{target_feature!r} not found in DataFrame columns.")

    # Compute correlations for the target feature
    corr = df.corr(method=method)[target_feature].drop(target_feature)

    # Decide sorting key
    sort_values = corr.abs() if abs_corr else corr
    corr_sorted = corr.loc[sort_values.sort_values(ascending=True).index]

    # Optionally keep only top_n by absolute magnitude
    if top_n is not None:
        corr_sorted = corr_sorted.iloc[-top_n:]

    # Create axis if needed
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, max(4, len(corr_sorted) * 0.25)))
        created_fig = True

    # Positions and values
    y_pos = list(range(len(corr_sorted)))
    x_vals = corr_sorted.values
    labels = corr_sorted.index

    # Draw stems (lines) and heads (dots)
    ax.hlines(y=y_pos, xmin=0, xmax=x_vals)
    ax.scatter(x_vals, y_pos)

    # ---- Padding + annotation ----
    if len(x_vals) > 0:
        x_min = float(min(min(x_vals), 0))
        x_max = float(max(max(x_vals), 0))
        x_range = x_max - x_min if x_max != x_min else 1.0

        # More generous padding than before, to give room for text
        left_pad = 0.45 * x_range
        right_pad = 0.45 * x_range
        ax.set_xlim(x_min - left_pad, x_max + right_pad)

        # Slightly bigger offset so text doesn't sit on the dots
        offset = 0.04 * x_range

        for y, x in zip(y_pos, x_vals):
            if x >= 0:
                text_x = x + offset
                ha = "left"
            else:
                text_x = x - offset
                ha = "right"

            ax.text(
                text_x,
                y,
                f"{x:.{decimals}f}",
                va="center",
                ha=ha,
                fontsize=label_fontsize,
            )
    # ------------------------------

    # Formatting
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel(f"{method.capitalize()} correlation with {target_feature}")
    # ax.set_title(f"Lollipop plot of correlations with {target_feature}")
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.margins(y=0.02)

    plt.tight_layout()
    if created_fig:
        plt.show()

    return ax


################################################################################
from scipy.signal import savgol_filter, find_peaks

# Define a function to smooth the data using a Savitzky-Golay filter
def smooth_data(data, window_length=5, polyorder=3):
    return savgol_filter(data, window_length, polyorder)


def find_end_slope2(data, height=None, plot=False, subject=None,smooth=False):

    if smooth:
        data = smooth_data(data)
    
    min_max_indc=[np.argmax(data), np.argmin(data)]
    last_peak_ind0=[i for i in min_max_indc if (i<len(data)-2)]# and (i>0)]
    if last_peak_ind0==[]:
        return 0
    
    last_peak_ind=np.max(last_peak_ind0)

    slope = (data[-1] - data[last_peak_ind]) / (len(data) - last_peak_ind-1)
    if plot:
        plt.figure()
        x_values = range(len(data))
        plt.plot(x_values, data, label="Data", color="blue",linestyle='-', marker='o')
        y_values_slope = [data[last_peak_ind] + slope * (x - last_peak_ind) for x in x_values]
        plt.ylim([-1.5,1.5])
        plt.plot(x_values, y_values_slope, label="Slope", color="red")
        plt.savefig(rootDir+"/workspace/results/slope_subjects_donna2/"+subject+'_smooth.png')
        
        plt.legend()
        plt.show()
#         plt.pause(0.01)

    return last_peak_ind, slope


