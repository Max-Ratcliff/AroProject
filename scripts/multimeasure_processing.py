# %%
import pandas as pd
from scipy.signal import savgol_filter
import numpy as np
import os
# import re
import math
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
import json
# %%
# def add_exp_labels(df, exp):
#     df["exp"] = [exp] * len(df)

#     exp_labels = {
#         "M4576_s2": "10mM Pulses #3",
#         "M4581_s1": "10mM Pulses #1",
#         "M4584_s1": "10mM Pulses #2",
#         "M6881_s2": "1mM Pulses #1",
#         "M6881_s5": "0.1mM Pulses #1",
#         "M6881_s6": "0.01mM Pulses #1",
#         # autotht
#         "M6813_s1": "0.01mM Cont. #1",
#         "M6813_s4": "0.01mM Cont. #2",
#         "M6605_s3": "10mM Pulses #1",
#         "M6605_s4": "10mM Pulses #2",
#         "M6605_s8": "10mM Pulses #3"
#     }

#     df["Exp_Label"] = df["exp"].map(exp_labels)
#     return df


def add_exp_labels(df, exp):
    df["exp"] = [exp] * len(df)

    exp_labels = {
        "Leticia_M4576_s2": "Leticia M4576 s2",
        "Pulses_M4581_s2": "Pulses M4581 s2",
        "Pulses_M4584_s1": "Pulses M4584 s1",
        # Add more descriptive names for your experiments here
    }

    df["Exp_Label"] = df["exp"].map(exp_labels)
    # If an experiment is not in the map, use its raw name as the label
    df["Exp_Label"].fillna(df["exp"], inplace=True)
    return df

# %%


def add_delta_intensity(df):
    df = df.copy()
    df["Delta_ThT"] = df["Intensity_ThT"] - df["Initial_Intensity_ThT"]
    return df

# %%


def plot_styling(df, column, cmap_name="tab20", num_samples=20):
    unique_vals = df[column].unique()
    unique_vals = list(unique_vals)[:num_samples]
    cmap = plt.get_cmap(cmap_name, len(unique_vals))
    color_dict = {val: cmap(i) for i, val in enumerate(unique_vals)}

    return color_dict


# %%


def add_germination_time(df, track_column, feature="DerivSavgol_Intensity", threshold=-8):
    spore_data = df.groupby(track_column)
    phase_df = []
    skip_ids = []
    suffix = track_column.split('_')[-1]  # Extracts 'PhC' or 'ThT'

    for spore_id, spore in spore_data:
        spore = spore.sort_values("Frame").copy()
        drop_frames = spore[spore[feature] < threshold]
        germ_frame = drop_frames["Frame"].min() if not drop_frames.empty else None

        # Create correctly named columns
        spore[f"Germination_Index_{suffix}"] = germ_frame

        if germ_frame is not None:
            # Compare every frame to the germination frame
            spore[f"Status_{suffix}"] = (spore["Frame"] >= germ_frame).astype(int)
        else:
            spore[f"Status_{suffix}"] = 0  # If it never germinates, status is always 0
            skip_ids.append(spore_id)

        phase_df.append(spore)

    if not phase_df:
        return pd.DataFrame()  # Return empty if no spores were processed

    df = pd.concat(phase_df, ignore_index=True)
    return df[~df[track_column].isin(skip_ids)]

# %%


def add_derivative_savgol(df, feature, track_column, window_length=7, poly_order=3):
    df = df.copy()
    deriv_col = f"DerivSavgol_{feature}"
    results = []

    spore_data = df.groupby([track_column])
    for spore_id, spore in spore_data:
        spore = spore.sort_values("Frame").copy()
        x = spore["Frame"].values
        y = spore[feature].values.astype(float)

        dt = np.median(np.diff(x))
        dy_dx = savgol_filter(y, window_length, poly_order, deriv=1, delta=dt)

        spore[deriv_col] = dy_dx
        results.append(spore)

    return pd.concat(results, ignore_index=True)

# %%


def extract_xy(df):
    pos_dict = {}
    spore_data = df.groupby(["Track_ID"])
    for track, spore in spore_data:
        avg_x = spore["X"].mean()
        avg_y = spore["Y"].mean()
        track_id = spore["Track_ID"].values[0]
        pos_dict[track_id] = (avg_x, avg_y)
    return pos_dict

# %%


def match_ids(dict1, dict2, max_distance=5):
    matches = []
    used_tht_ids = set()

    for id1, (x1, y1) in dict1.items():
        closest_id = None
        min_dist = float('inf')

        for id2, (x2, y2) in dict2.items():
            if id2 in used_tht_ids:
                continue  # skip if already matched

            dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)
            if dist < min_dist and dist <= max_distance:
                min_dist = dist
                closest_id = id2

        if closest_id is not None:
            matches.append((id1, closest_id))
            used_tht_ids.add(closest_id)

    print(f"matched {len(matches)} tracks...")
    return matches


def merge_matched_tracks(phc_df, tht_df, matched_pairs, phc_suffix="PhC", tht_suffix="ThT"):
    merged_list = []
    for phc_id, tht_id in matched_pairs:
        phc_track = phc_df[phc_df["Track_ID"] == phc_id].copy()
        tht_track = tht_df[tht_df["Track_ID"] == tht_id].copy()

        phc_track = phc_track.add_suffix(f"_{phc_suffix}")
        tht_track = tht_track.add_suffix(f"_{tht_suffix}")

        phc_track = phc_track.rename(columns={f"Frame_{phc_suffix}": "Frame"})
        tht_track = tht_track.rename(columns={f"Frame_{tht_suffix}": "Frame"})

        merged = pd.merge(phc_track, tht_track, on="Frame", how="right")
        # merged["MATCHED_PAIR"] = f"{phc_id}_{tht_id}"
        merged_list.append(merged)

    if not merged_list:
        return pd.DataFrame()

    return pd.concat(merged_list, ignore_index=True)


# %%


def calculate_interval_germ_exposures(initial_min, between_min, minutes_between_frames, exp_length):
    i = initial_min/minutes_between_frames - 1
    exposures = [i]
    while i < exp_length:
        i += between_min/minutes_between_frames
        exposures.append(i)
    if exposures[-1] >= exp_length:
        exposures = exposures[:-1]
    return exposures


def add_germ_exposures_todf(df, germ_exposures_list, concentration):
    is_exposure = df["Frame"].isin(germ_exposures_list)
    df["Germinant"] = is_exposure.astype(int) * concentration
    return df

# %%


def plot_intensity(df, track_column, color_dict, feature="Intensity_ThT", alpha=1, linestyle="-", num_samples=None):
    germ_col = f"Status_{track_column.split('_')[-1]}"
    spore_data = df.groupby(track_column)
    if num_samples is not None:
        samples = 0
    for track_id, data in spore_data:
        # exp = data["exp"].iloc[0]  # Add this
        color = color_dict.get(track_id, "gray")

        dormant_data = data[data[germ_col] == 0]
        germinated_data = data[data[germ_col] == 1]

        sns.lineplot(x="Frame", y=feature, data=dormant_data,
                     linewidth=2, color=color, alpha=alpha, label=feature.replace("_", " "), linestyle=linestyle)
        # sns.lineplot(x = "FRAME", y = feature, data = data, alpha = alpha, color = color)
        if not germinated_data.empty:
            sns.lineplot(x="Frame", y=feature, data=germinated_data,
                         linewidth=6, color=color, alpha=alpha, linestyle=linestyle)
        if samples > num_samples - 1:
            return
        samples += 1


# %%


def manually_filter_spores(base, df, track_column, feature, display_time=1, filter=1):

    # Grab experiment name safely
    exp = df["exp"].iloc[0] if "exp" in df.columns else "unknown"

    # Set up ignore file path
    ignore_file = os.path.join(base, f"{exp}_FilteredPhCTracks.csv")

    # Load ignored spores
    if os.path.exists(ignore_file):
        ignored_spores = set(map(str, pd.read_csv(ignore_file)[track_column]))
    else:
        ignored_spores = set()

    spores_to_keep = []
    spores_to_ignore = []
    if filter == 0:
        print("skipping filtering...")
        print(f'keeping {df[track_column][~df[track_column].isin(ignored_spores)].nunique()} spores...')
        return df[~df[track_column].isin(ignored_spores)]
    spore_data = df.groupby(track_column)
    for spore_id, data in spore_data:

        if spore_id in ignored_spores:
            print("found")
            continue

        # Plot
        fig, ax = plt.subplots(figsize=(3, 1))
        sns.lineplot(data=data, x="Frame", y=feature, ax=ax)
        ax.axvline(data["Germination_Index_" + track_column.split('_')[-1]].values[0], color='red', linestyle='--')
        plt.show(block=False)
        plt.pause(display_time)
        plt.close(fig)

        keep = input(f"Keep spore {spore_id}? (y/n/a): ").strip().lower()
        if keep not in ["a", "y", "n"]:
            print("Invalid input... showing again...")
            keep = "a"

        if keep == "a":
            fig, ax = plt.subplots(figsize=(3, 1))
            sns.lineplot(data=data, x="Frame", y=feature, ax=ax)
            ax.axvline(data["Germination_Index_PhC"].values[0], color='red', linestyle='--')
            plt.show(block=False)
            plt.pause(display_time)
            plt.close(fig)
            keep = input(f"Keep spore {spore_id}? (y/n): ").strip().lower()

        if keep == "y":
            spores_to_keep.append(spore_id)
        elif keep == "n":
            spores_to_ignore.append(spore_id)

    # Save ignored spores
    if spores_to_ignore:
        new_ignored = pd.DataFrame(spores_to_ignore, columns=[track_column])
        if os.path.exists(ignore_file):
            existing = pd.read_csv(ignore_file)
            updated = pd.concat([existing, new_ignored], ignore_index=True).drop_duplicates()
        else:
            updated = new_ignored
        updated.to_csv(ignore_file, index=False)

    print(f"Keeping {len(spores_to_keep)} spores...")

    # Return filtered DataFrame
    return df[df["Track_PhC"].isin(spores_to_keep)]


# %%
def add_rolling_derivative(df, feature, track_column, window=5):
    df = df.copy()
    deriv_col = f"Derivative_{feature}"
    results = []

    spore_data = df.groupby(["Exp_Label", track_column])
    for (exp_label, spore_id), spore in spore_data:
        spore = spore.sort_values("Frame").copy()
        rolling_diff = spore[feature].diff().rolling(window=window, min_periods=1).mean()
        spore[deriv_col] = rolling_diff.fillna(0)
        results.append(spore)

    return pd.concat(results, ignore_index=True)


def add_rolling_mean(df, feature, track_column, window=5):
    df = df.copy()
    smoothed_col = f"Rolling_{feature}"
    result = []

    spore_data = df.groupby(["Exp_Label", track_column])
    for spore_id, spore in spore_data:
        spore = spore.sort_values("Frame").copy()
        spore[smoothed_col] = spore[feature].rolling(window=window, min_periods=1).mean()
        result.append(spore)

    return pd.concat(result, ignore_index=True)

# %%


def add_initial_feature(df, feature, track_column, num_init=11):
    df = df.copy()
    df = df.sort_values("Frame")

    init_col_name = f"Initial_{feature}"
    column_values = []

    spore_data = df.groupby(["Exp_Label", track_column])
    for spore_id, spore in spore_data:
        feature_vals = spore[feature].values.astype(float)
        initial_feature = np.mean(feature_vals[:num_init])

        spore[init_col_name] = [initial_feature] * len(spore)
        column_values.append(spore)

    df_return = pd.concat(column_values, ignore_index=True)
    return df_return


def add_cumsum_feature(df, feature, track_column):
    df = df.sort_values([track_column, "Frame"]).copy()
    cumsum_col = f"CumuSum_{feature}"
    df[cumsum_col] = df.groupby(track_column)[feature].cumsum()
    return df

# %%


def plot_phase(df, track_column, feature="DerivSavgol_Intensity_PhC", frame_col="Frame", intensity_col="Intensity_PhC"):
    for spore_id, spore_data in df.groupby(track_column):
        # plt.axhline(-8, color = "lightgrey", label = "Threshold")
        # germinated_df = spore_data[spore_data["Status_PhC"] == 1]
        # dormant_df = spore_data[spore_data["Status_PhC"] == 0]
        sns.lineplot(x=frame_col, y=intensity_col, data=spore_data, color="tab:blue", linewidth=3, label="Phase Intensity", alpha=1)

        sns.lineplot(x=frame_col, y=feature, data=spore_data, color="tab:orange", linewidth=3, label="Savitzky-Golay Derivative")
        # sns.lineplot(x=frame_col, y=feature, data=spore_data, color="tab:orange", linewidth=1, alpha=0.5)
        # sns.lineplot(x=frame_col, y=feature, data=germinated_df, color="orange", linewidth=5)
        plt.axvline(spore_data["Germination_Index_PhC"].values[0], label="Threshold Met", linewidth=3, color="darkgrey", alpha=0.7)

        # sns.lineplot(x=frame_col, y=intensity_col, data=germinated_df, color="tab:blue", linewidth=5, label="PhC: phase-dark")
        plt.xlabel("Frame", fontsize=14)
        plt.ylabel("")
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.legend(fontsize=12)
        plt.show()


# %%


def add_exposures_during_dormancy(df, track_column, germ_exposures_list):
    """Add count of germinant exposures each spore experienced while dormant"""
    results = []

    for spore_id, spore_data in df.groupby(track_column):
        spore_data = spore_data.sort_values("Frame").copy()
        germ_frame = spore_data["Germination_Index_PhC"].iloc[0]

        if pd.isna(germ_frame):
            # Never germinated - count all exposures
            exposures_during_dormancy = len([exp for exp in germ_exposures_list if exp <= spore_data["Frame"].max()])
        else:
            # Count exposures before germination
            exposures_during_dormancy = len([exp for exp in germ_exposures_list if exp < germ_frame])

        spore_data["Exposures_Dormancy"] = exposures_during_dormancy
        results.append(spore_data)

    return pd.concat(results, ignore_index=True)

# %%


# --- Configuration ---
config_path = Path(__file__).resolve().parents[1] / 'config.json'
with open(config_path, 'r') as f:
    config = json.load(f)

DATA_ROOT = config["DATA_ROOT"]
EXPERIMENT = config["CURRENT_EXPERIMENT"]

if not all([DATA_ROOT, EXPERIMENT]):
    raise ValueError("DATA_ROOT or CURRENT_EXPERIMENT not set in your config.json file.")

# This variable is for labeling the output file, you can change it as needed
EXP_LABEL = f"{EXPERIMENT}_Analysis"

# --- Static variables for this script's specific purpose ---
# This script's job is to merge PhC and a fluorescent channel.
MICR_PHC = "PhC"
MICR_FLUOR = "ThT"

BASE = Path(DATA_ROOT)

# --- Path Setup ---
# CORRECTED: Added "data" and "processed" to the path to match your structure
fiji_base = BASE / "data" / "processed" / f"{EXPERIMENT}_Fiji"
data_folder = fiji_base / "Processed_Data"

# Create the output directory if it doesn't exist to prevent errors on save.
data_folder.mkdir(parents=True, exist_ok=True)

# --- Data Loading and Initial Processing ---
# Initialize dataframes as None. We will check later if they were successfully loaded.
phc_df = None
tht_df = None

print(f"Searching for results files in: {fiji_base}")
try:
    for csv_filename in os.listdir(fiji_base):
        if "_Results.csv" not in csv_filename:
            print(f"Skipping non-results file: {csv_filename}")
            continue  # Skip any file that isn't a results file

        # Full path to the current CSV file
        input_csv_path = os.path.join(fiji_base, csv_filename)

        if f"_{MICR_PHC}_" in csv_filename:
            print("-> Found 'PHC' data, loading...")
            phc_df = pd.read_csv(input_csv_path)
        elif f"_{MICR_FLUOR}_" in csv_filename:
            print("-> Found 'ThT' data, loading...")
            tht_df = pd.read_csv(input_csv_path)

except FileNotFoundError:
    print("---")
    print("ERROR: The directory was not found.")
    print("Please verify this path exists: {fiji_base}")
    phc_df, tht_df = None, None  # Ensure dataframes are None so merge step is skipped


# --- Conditional Merging and Analysis ---
# This block only runs if BOTH phc_df and tht_df were successfully loaded in the step above.
if phc_df is not None and tht_df is not None:
    print("\nSUCCESS: Both PhC and ThT data found. Merging and analyzing...")

    # Match spores between the two channels based on XY coordinates
    phc_xy = extract_xy(phc_df)
    tht_xy = extract_xy(tht_df)
    matches = match_ids(phc_xy, tht_xy)

    # Merge the dataframes into a single master dataframe 'df'
    df = merge_matched_tracks(phc_df, tht_df, matches, phc_suffix="PhC", tht_suffix="ThT")
    df = add_exp_labels(df, EXPERIMENT)

    # --- Feature Engineering and Analysis on Merged Data ---
    # Process PhC data
    TRACK_COL_MERGED = "Track_ID_PhC"
    df = add_derivative_savgol(df, "Intensity_PhC", TRACK_COL_MERGED, window_length=7, poly_order=3)
    print("columns:", df.columns)
    df = add_germination_time(df, TRACK_COL_MERGED, feature="DerivSavgol_Intensity_PhC")

    # Check if germination detection worked
    germ_count = df["Germination_Index_PhC"].notna().sum()
    print(f"Detected germination in {germ_count} spore tracks")

    if germ_count == 0:
        print("WARNING: No germination events detected. Check threshold and feature values.")
        print(f"DerivSavgol_Intensity_PhC range: {df['DerivSavgol_Intensity_PhC'].min()} to {df['DerivSavgol_Intensity_PhC'].max()}")

    # Add germinant exposures (you'll need to define these based on your experiment)
    germ_exposures_list = calculate_interval_germ_exposures(
        initial_min=10, between_min=30, minutes_between_frames=1, exp_length=289
    )

    # Add exposures during dormancy column
    df = add_exposures_during_dormancy(df, TRACK_COL_MERGED, germ_exposures_list)

    df = add_initial_feature(df, "Intensity_PhC", TRACK_COL_MERGED)

    # Process ThT data
    feature_tht = f"Intensity_{MICR_FLUOR}"
    df = add_rolling_mean(df, feature_tht, TRACK_COL_MERGED)
    df = add_rolling_derivative(df, feature_tht, TRACK_COL_MERGED)
    df = add_initial_feature(df, feature_tht, TRACK_COL_MERGED)
    df = add_cumsum_feature(df, feature_tht, TRACK_COL_MERGED)
    df = add_delta_intensity(df)

    # Plot the results
    # plot_phase(df, TRACK_COL_MERGED)

    # Save the final, merged, and processed data
    final_output_path = os.path.join(data_folder, f"{EXP_LABEL}_Matched_And_Processed_Data.csv")
    df.to_csv(final_output_path, index=False)
    print(f"\nANALYSIS COMPLETE. Final merged data saved to:\n{final_output_path}")

else:
    print("---")
    print("SKIPPING MERGE STEP: Could not load both PhC and ThT data.")
