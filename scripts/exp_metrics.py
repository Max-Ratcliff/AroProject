# %%
# Libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
# import ast
# import math
import seaborn as sns
# from sklearn.preprocessing import MinMaxScaler
from matplotlib.ticker import FuncFormatter
# from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
from pathlib import Path
from dotenv import load_dotenv
# from sklearn.cross_decomposition import PLSRegression
# import matplotlib.colors as mcolors

# import matplotlib.cm as cm
# from sklearn.decomposition import PCA
# from sklearn.preprocessing import StandardScaler

# %%


def get_bin_width(data, method):
    if method == "Sturges":
        k = int(np.ceil(np.log2(len(data)) + 1))
        bin_width = (np.max(data) - np.min(data)) / k

    if method == "Scotts":
        bin_width = 3.5 * np.std(data) / (len(data) ** (1/3))

    if method == "Freedman":
        q25, q75 = np.percentile(data, [25, 75])
        iqr = q75 - q25
        bin_width = 2 * iqr / (len(data) ** (1/3))

    bins = int(np.ceil((data.max() - data.min()) / bin_width))
    return bins

# %%


def plot_styling(df, track_id_col, cmap_name="tab20", num_samples=None):
    # Grab the unique combinations of 'exp' and track_id
    df_unique = df[["exp", track_id_col]].drop_duplicates()

    # Optionally limit the number of samples
    if num_samples is not None:
        df_unique = df_unique.head(num_samples)

    # Turn into list of tuples for consistent color assignment
    unique_ids = list(df_unique.itertuples(index=False, name=None))

    # Build color dictionary from colormap
    cmap = plt.get_cmap(cmap_name, len(unique_ids))
    color_dict = {uid: cmap(i) for i, uid in enumerate(unique_ids)}

    return color_dict

# %%


def histogram_germination_metrics(plots_folder: str, df, track_col, germinant_exposures=None, output=None, title=None) -> None:
    plt.figure(figsize=(6, 3))
    palette = sns.color_palette("colorblind", n_colors=2)

    spore_germ_frames = df.groupby(track_col)["Germination_Index_PhC"].first()
    # 2. total spores
    total_spores = spore_germ_frames.shape[0]

    # 3. spores that actually germinated
    germinated_spores = spore_germ_frames.dropna()

    # 5. list of germination frames
    germination_frames_list = germinated_spores.tolist()
    last_germ_time = 289  # you manually set it

    frames_shown: int = 289
    bin_size = get_bin_width(df[track_col].unique(), "Sturges")


# Auto-infer germinant exposures if not given
    if not germinant_exposures:
        germinant_exposures = df.loc[df["Germination_Index_PhC"] != 0, 'Frame'].unique()
        # print(f"found germinant exposures {germinant_exposures}")

    # -- Germination event counting
    frame_counts = Counter(germination_frames_list)
    sorted_frame_counts = sorted(frame_counts.items())

    germination_events: list[int] = []
    percent_germinated_at_t: list[int] = []

    for frame_number, count in sorted_frame_counts:
        if output == 1:
            germination_events.append(count)

    # -- Dormant percentage over time
    # total_spores = sum(count for frame, count in sorted_frame_counts)
    spores_count = total_spores  # (do NOT overwrite total_spores!)

    frames = [0]
    percents = [100]
    percent_plot = [100]

    frame_dict = dict(sorted_frame_counts)
    total_percent = 100

    for frame_number in range(1, frames_shown + 1):
        if frame_number in frame_dict:
            count = frame_dict[frame_number]
            spores_count -= count
            total_percent = spores_count / total_spores * 100
            t_percent = count / total_spores * 100
            percent_plot.append(total_percent)
            percent_germinated_at_t.append(t_percent)
        percents.append(total_percent)
        frames.append(frame_number)

    fig, ax2 = plt.subplots(figsize=(6, 4))  # <-- make ax2 first now (flip)
    palette = sns.color_palette("colorblind", n_colors=2)

    # Plot dormant % first, on LEFT (primary)
    sns.lineplot(x=frames, y=percents, ax=ax2, linewidth=4, label="Dormant Percentage", color=palette[0])
    ax2.set_ylabel("Population Dormancy", fontsize=16)
    ax2.set_yticks([0, 25, 50, 75, 100])
    ax2.set_ylim([-10, 110])  # **use set_ylim, not set_ybound**
    ax2.tick_params(axis='y', labelsize=12)
    ax2.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{int(y)}%'))
    ax2.set_xlim([0, frames_shown])

    # Create secondary axis for RIGHT (for germination count bars)
    ax1 = ax2.twinx()

    sns.histplot(germination_frames_list, bins=bin_size, label="Germination Events", ax=ax1, color=palette[1], alpha=0.5)
    ax1.set_ylabel("Germination Count", fontsize=16)
    ax1.set_xlabel("Frame", fontsize=16)
    ax1.set_xlim([0, frames_shown])
    # -- Legends

    # -- Plot germinant exposures if they exist
    if germinant_exposures.all() is not None:
        for exposure in germinant_exposures:
            ax2.axvline(x=exposure, color='silver', linewidth=1, linestyle="dashed", label='Germinant Exposure' if exposure == germinant_exposures[0] else "")

    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(handles=handles1 + handles2, labels=labels1 + labels2, loc='best', fontsize=12)

    ax1.set_xlim([0, last_germ_time])
    ax2.set_xlim([0, last_germ_time])
    if title is None:
        title = df["Exp_Label"].values[0]
    plt.title(title, fontsize=22)
    if title == "0.01mM Pulses Exp. 6":
        ax2.legend(handles=handles1 + handles2, labels=labels1 + labels2, loc='lower left', fontsize=12)

    plt.tight_layout()

    # If you want to save it:
    # plt.savefig(f"{plots_folder}hist_germinationmetrics_{bin_size}bins.jpg")
    # plt.clf()

    # -- Stats
    if output:
        print(f"{total_spores} total spores...")

        print("\nStatistics by Germination Frames:")
        print(f"Frames that events occur: {len(germination_events)}")
        print(f"Frames: {str(list(frame_dict.keys())).replace(',', ' &')}")
        print(f"Germination Events: {str(germination_events).replace(',', ' &')}")
        print(f"Total Percentage: [{' & '.join([f'{elem:.0f}' for elem in percent_plot[1:]])}]")
        print(f"Percentage at frame: [{' & '.join([f'{elem:.0f}' for elem in percent_germinated_at_t])}]")

    return None


# %%

def lineplot_feature_trends(df, feature: str, track_col, max_frame=289, title=None):
    fig, ax = plt.subplots(figsize=(6, 4))

    exposure_groups = df.groupby("Exposures_Dormancy")
    colors = sns.color_palette("tab20", n_colors=exposure_groups.ngroups)

    for idx, (exposures, exposure_group) in enumerate(exposure_groups):
        color = colors[idx]

        if exposures > -1:
            dormant = exposure_group[exposure_group["Status_PhC"] == 0]
            germinated = exposure_group[exposure_group["Status_PhC"] == 1]

            # Plot dormant phase (thin line)
            sns.lineplot(
                x="Frame",
                y=feature,
                data=dormant,
                ax=ax,
                color=color,
                linewidth=3,
                label=None,  # no label for thin line,
                alpha=0.5,
                errorbar='se',
            )

            # Plot germinated phase (thick line)
            sns.lineplot(
                x="Frame",
                y=feature,
                data=germinated,
                ax=ax,
                color=color,
                linewidth=7,
                label=f"{int(exposures)}",
                alpha=0.8,
                errorbar='se'
            )
        else:
            # For dormant spores that never germinated
            sns.lineplot(
                x="Frame",
                y=feature,
                data=exposure_group,
                ax=ax,
                color=color,
                linewidth=2,
                label="",
                alpha=0.5,
                errorbar='se'
            )

    ax.set_xlim([0, max_frame])
    ax.set_xlabel("Frame", fontsize=18)
    ax.set_ylabel(f'{feature.replace("_", " " ).replace("ThT", ""). replace("CH1", "").title()}', fontsize=18)
    # ax.legend(title="Exposures during Dormancy")
    ax.set_title(title, fontsize=20)
    plt.tight_layout()


# %%


def scatterplot_initial_autotht(df, track_col):
    df_initial = df[df["FRAME"] == 0]
    sns.scatterplot(data=df_initial, x="INIT_ROLLAVG_MEAN_INTENSITY_CH1_ThT_mCHE", y="INIT_ROLLAVG_MEAN_INTENSITY_CH1_ThT_CFP", hue="GERMINATION_FRAME", palette="pastel")

# %%


def lineplot_feature(df, feature, color_dict, germ_exp=None, track_column="TRACK_ID_ThT", num_samples=None, alpha=1):

    spore_data = df.groupby(track_column)
    spore_count = 0
    for track_id, data in spore_data:
        exp = data["exp"].iloc[0]  # Add this
        color = color_dict.get((exp, track_id), "gray")
        if "Germination_Index_PhC" in data.columns:
            germ_frame = data["Germination_Index_PhC"].iloc[0]
            dormant_data = data[data["Frame"] < germ_frame]
            germinated_data = data[data["Frame"] >= germ_frame]
        else:
            dormant_data = data
            germinated_data = pd.DataFrame()

        sns.lineplot(x="Frame", y=feature, data=dormant_data,
                     linewidth=2, color=color, alpha=alpha)
        sns.lineplot(x="Frame", y=feature, data=data, alpha=0.5, color=color)
        if not germinated_data.empty:
            sns.lineplot(x="Frame", y=feature, data=germinated_data,
                         linewidth=6, color=color, alpha=alpha)

        spore_count += 1
        if num_samples:
            if spore_count >= num_samples:
                break
    if germ_exp:
        for germ in germ_exp:
            plt.axvline(germ, alpha=0.4, color="lightgrey")
    plt.xlabel("Frame")
    plt.ylabel(feature.replace("_", " "), fontsize=16)
    # plt.ylim([-10, 200])


# %% [markdown]
# ---

# %%
load_dotenv()
DATA_ROOT = os.getenv("DATA_ROOT")
if DATA_ROOT is None:
    raise ValueError("DATA_ROOT environment variable is not set.")

# --- Analysis Constants ---
FEATURE_TO_ANALYZE = "Intensity_ThT"
TRACKING_COLUMN = "Track_ID_ThT"

# SETUP PATHS USING PATHLIB
PROJECT_ROOT = Path(DATA_ROOT)
PROCESSED_DATA_DIR = PROJECT_ROOT / "data" / "processed"
PLOTS_OUTPUT_DIR = PROJECT_ROOT / "reports" / "figures"

# Create the output directory for plots if it doesn't exist.
PLOTS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"Project Root: {PROJECT_ROOT}")
print(f"Data will be read from: {PROCESSED_DATA_DIR}")
print(f"Plots will be saved to: {PLOTS_OUTPUT_DIR}")

# rglob() recursively finds all files ending with "_Processed_Data.csv".
# pattern: processed/{exp_name}_Fiji/Processed_Data/{exp_name}_{phase}_Processed_Data.csv
data_files = list(PROCESSED_DATA_DIR.rglob("*/Processed_Data/*_Processed_Data.csv"))

if not data_files:
    print(f"\nWarning: No data files found in {PROCESSED_DATA_DIR}")
else:
    print(f"\nFound {len(data_files)} data files to process.")


for file_path in data_files:
    print(f"\n--- Processing: {file_path.name} ---")

    df = pd.read_csv(file_path)

    print(f"Columns found: {df.columns.to_list()}")

    experiment_name = file_path.name.replace("_Results.csv", "")  # grab experiment name

    plot_title = f"{experiment_name} ({FEATURE_TO_ANALYZE})"

    hist_output_path = PLOTS_OUTPUT_DIR / f"{experiment_name}_germination_histogram.png"
    trends_output_path = PLOTS_OUTPUT_DIR / f"{experiment_name}_feature_trends.png"

    # Plot 1: Histogram of Germination Metrics
    print("Generating histogram...")
    histogram_germination_metrics(PLOTS_OUTPUT_DIR, df, TRACKING_COLUMN, title=plot_title)
    plt.savefig(hist_output_path, bbox_inches='tight')
    plt.close()   # Close plot to free memory
    print(f"Saved plot: {hist_output_path}")

    # Plot 2: Lineplot of Feature Trends
    print("Generating lineplot...")
    lineplot_feature_trends(df, FEATURE_TO_ANALYZE, TRACKING_COLUMN, title=plot_title)
    plt.savefig(trends_output_path, bbox_inches='tight')
    plt.close()   # Close plot to free memory
    print(f"Saved plot: {trends_output_path}")

print("\nProcessing complete")
