import pandas as pd
# import numpy as np
from pathlib import Path
from skimage import io, measure
from natsort import natsorted
from tqdm import tqdm


def process_timelapse_with_mask(data_dir: Path, mask_file: Path, output_csv: Path):
    """
    Applies a single mask to a timelapse series of images, extracts spore properties,
    and saves the data in a wide format compatible with the processing notebook.

    Args:
        data_dir (Path): Directory containing the TIFF image sequence.
        mask_file (Path): Path to the Cellpose-generated '_cp_masks.png' file.
        output_csv (Path): Path to save the resulting CSV file.
    """
    # 1. Load the single mask file.
    try:
        mask_image = io.imread(mask_file)
        print(f"Successfully loaded mask file from: {mask_file}")
    except FileNotFoundError:
        print(f"Error: Mask file not found at {mask_file}")
        return

    # Find all image files in the data directory and sort them numerically.
    image_files = natsorted(
        [f for f in data_dir.glob("*.tif") if "_masks" not in f.name]
    )
    if not image_files:
        print(f"Error: No .tif files found in {data_dir}")
        return

    print(f"Found {len(image_files)} images to process in {data_dir}")

    all_frames_data = []

    # 2. Iterate through each image file in the timelapse.
    print("Processing frames...")
    for frame_index, image_path in enumerate(tqdm(image_files, desc="Analyzing Frames")):
        current_frame_image = io.imread(image_path)

        # 3. Apply the mask to the current frame to measure properties.
        properties = measure.regionprops_table(
            mask_image,
            intensity_image=current_frame_image,
            properties=('label', 'mean_intensity', 'centroid', 'area', 'perimeter', 'major_axis_length', 'minor_axis_length')
        )

        frame_df = pd.DataFrame(properties)

        # Format columns to match your `extract_features` function
        frame_dict = {'Slice': frame_index + 1}
        for _, row in frame_df.iterrows():
            spore_id = row['label']
            frame_dict[f'Mean{spore_id}'] = row['mean_intensity']
            frame_dict[f'X{spore_id}'] = row['centroid-1']  # centroid-1 is X
            frame_dict[f'Y{spore_id}'] = row['centroid-0']  # centroid-0 is Y
            frame_dict[f'Area{spore_id}'] = row['area']
            frame_dict[f'Perim.{spore_id}'] = row['perimeter']
            frame_dict[f'Major{spore_id}'] = row['major_axis_length']
            frame_dict[f'Minor{spore_id}'] = row['minor_axis_length']

        all_frames_data.append(frame_dict)

    # 4. Create the final wide-format DataFrame and save to CSV.
    final_df = pd.DataFrame(all_frames_data)
    final_df.to_csv(output_csv, index=False)
    print(f"\nProcessing complete. Data saved to {output_csv}")


# --- Configuration ---
BASE_DIR = Path("/Users/mratcliff/Desktop/AroProject/")
EXP_FOLDER = "Leticia_M4576_s2"
MICROSCOPY_TYPE = "PH"

# The source directory for the original images
DATA_DIR = BASE_DIR / EXP_FOLDER / "M4576_s2_PH"

# The output directory for the results CSV
OUTPUT_DIR = BASE_DIR / f"{EXP_FOLDER}_Fiji"
OUTPUT_DIR.mkdir(exist_ok=True)

# --- MODIFIED LINE ---
# Directly point to the mask file generated in your '/test/' folder.
MASK_FILE = Path("/Users/mratcliff/Desktop/AroProject/test/M4576_s2_PH_stabilized_0000_cp_masks.png")

# Construct the output CSV filename as required by your processing script.
# Using "PhC" for compatibility with your notebook.
OUTPUT_CSV = OUTPUT_DIR / f"{EXP_FOLDER}_PhC_Results.csv"


# --- Run Processing ---
process_timelapse_with_mask(DATA_DIR, MASK_FILE, OUTPUT_CSV)
