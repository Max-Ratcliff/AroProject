import pandas as pd
# import numpy as np
from pathlib import Path
from skimage import io, measure
from natsort import natsorted
from tqdm import tqdm
from dotenv import load_dotenv
import os


def process_timelapse_with_mask(data_dir: Path, mask_file: Path, output_csv: Path):
    """
    Applies a single mask to a timelapse, extracts spore properties,
    and saves the data in a robust, long format. This version is functionally
    equivalent to the original wide-format script but produces a more usable output.

    Args:
        data_dir (Path): Directory containing the TIFF image sequence.
        mask_file (Path): Path to the Cellpose-generated '_cp_masks.png' file.
        output_csv (Path): Path to save the resulting long-format CSV file.
    """

    try:
        mask_image = io.imread(mask_file)
        print(f"Successfully loaded mask file: {mask_file.name}")
    except FileNotFoundError:
        print(f"Error: Mask file not found at {mask_file}")
        return

    image_files = natsorted(
        [f for f in data_dir.glob("*.tif") if "_masks" not in f.name]
    )

    if not image_files:
        print(f"Error: No .tif files found in {data_dir}")
        return
    print(f"Found {len(image_files)} images to process in {data_dir}")

    all_frames_data = []

    # Iterate through each image file, showing progress with tqdm
    print("Extracting measurements from frames...")
    for frame_index, image_path in enumerate(tqdm(image_files, desc="Analyzing Frames")):
        current_frame_image = io.imread(image_path)

        # The core measurement step is identical to your original script
        properties = measure.regionprops_table(
            mask_image,
            intensity_image=current_frame_image,
            properties=('label', 'mean_intensity', 'centroid', 'area', 'perimeter', 'major_axis_length', 'minor_axis_length')
        )

        # if no spores are found, an empty DataFrame is created.
        frame_df = pd.DataFrame(properties)

        frame_df['Frame'] = frame_index + 1
        all_frames_data.append(frame_df)

    # Combine data from all frames into a single long-format DataFrame
    final_df = pd.concat(all_frames_data, ignore_index=True)

    # Rename columns for clarity and to remove special characters (like '-')
    final_df.rename(columns={
        'label': 'Track_ID',
        'mean_intensity': 'Intensity',
        'centroid-1': 'X',
        'centroid-0': 'Y',
        'area': 'Area',
        'perimeter': 'Perimeter',
        'major_axis_length': 'Major',
        'minor_axis_length': 'Minor'
    }, inplace=True)

    # Ensure the directory to save the file exists
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    final_df.to_csv(output_csv, index=False)
    print(f"\nProcessing complete. Long-format data saved to {output_csv}")


# --- Configuration ---
load_dotenv()  # Load environment variables from .env file if needed
# Load all configuration from the .env file
DATA_ROOT = os.getenv("DATA_ROOT")
EXPERIMENT = os.getenv("EXPERIMENT_NAME")
RAW_SUBFOLDER = os.getenv("RAW_SUBFOLDER")
MICROSCOPY_TYPE = os.getenv("MICROSCOPY_TYPE")
MASK_FILENAME = os.getenv("MASK_FILENAME")

# Check that variables were loaded successfully
if not all([DATA_ROOT, EXPERIMENT, RAW_SUBFOLDER, MICROSCOPY_TYPE, MASK_FILENAME]):
    raise ValueError("One or more required variables are not set in your .env file.")

DATA_ROOT = Path(DATA_ROOT)

# The source directory for the original images
DATA_DIR = DATA_ROOT / "data" / "raw" / EXPERIMENT / RAW_SUBFOLDER

# The output directory for the results CSV
OUTPUT_DIR = DATA_ROOT / "data/processed" / f"{EXPERIMENT}_Fiji"
OUTPUT_DIR.mkdir(exist_ok=True)

MASK_FILE = DATA_ROOT / "data/processed/cellpose" / f"{EXPERIMENT}_{MICROSCOPY_TYPE}" / MASK_FILENAME

# Construct the output CSV filename as required by your processing script.
# Using "PhC" for compatibility with your notebook.
OUTPUT_CSV = OUTPUT_DIR / f"{EXPERIMENT}_{MICROSCOPY_TYPE}_Results.csv"


# --- Run Processing ---
process_timelapse_with_mask(DATA_DIR, MASK_FILE, OUTPUT_CSV)
