import argparse
import json
from pathlib import Path

# import numpy as np
import pandas as pd
from cellpose import core, io, models
from skimage import measure


def process_frame(image_path: Path, config: dict):
    """
    Processes a single image file: runs Cellpose, measures properties, and appends to a CSV.
    """
    print(f"\n--- Processing new frame: {image_path.name} ---")

    # --- 1. Load Configuration ---
    DATA_ROOT = Path(config["DATA_ROOT"])
    EXPERIMENT = config["CURRENT_EXPERIMENT"]
    MICROSCOPY_TYPE = config["CURRENT_MICROSCOPY_TYPE"]

    # Define output paths based on config
    cellpose_output_dir = DATA_ROOT / f"data/processed/cellpose/{EXPERIMENT}_{MICROSCOPY_TYPE}"
    fiji_output_dir = DATA_ROOT / "data/processed" / f"{EXPERIMENT}_Fiji"
    results_csv_path = fiji_output_dir / f"{EXPERIMENT}_{MICROSCOPY_TYPE}_Results.csv"

    # Ensure directories exist
    cellpose_output_dir.mkdir(parents=True, exist_ok=True)
    fiji_output_dir.mkdir(parents=True, exist_ok=True)

    # --- 2. Run Cellpose for the single image ---
    print("Running Cellpose...")
    model = models.CellposeModel(gpu=core.use_gpu())
    img = io.imread(image_path)

    # model.eval returns masks, flows, styles
    masks, _, _ = model.eval([img], diameter=None, channels=[0, 0])
    mask_image = masks[0]  # Get the mask for our single image

    # Save the individual mask for this frame for debugging/review
    mask_output_path = cellpose_output_dir / f"{image_path.stem}_cp_masks.png"
    io.imsave(mask_output_path, mask_image)
    print(f"Saved individual mask to: {mask_output_path.name}")

    # --- 3. Measure Properties using the new mask ---
    print("Measuring properties...")
    properties = measure.regionprops_table(
        mask_image,
        intensity_image=img,
        properties=('label', 'mean_intensity', 'centroid', 'area', 'perimeter', 'major_axis_length', 'minor_axis_length')
    )
    frame_df = pd.DataFrame(properties)

    # Add Frame number from filename (e.g., '..._0001.tif' -> 2)
    try:
        frame_number = int(image_path.stem.split('_')[-1]) + 1  # Assuming 0-indexed filenames
        frame_df['Frame'] = frame_number
    except (ValueError, IndexError):
        print(f"Warning: Could not parse frame number from {image_path.name}. Skipping append.")
        return

    # Rename columns to match your existing schema
    frame_df.rename(columns={
        'label': 'Track_ID',
        'mean_intensity': 'Intensity',
        'centroid-1': 'X',
        'centroid-0': 'Y',
        'area': 'Area',
        'perimeter': 'Perimeter',
        'major_axis_length': 'Major',
        'minor_axis_length': 'Minor'
    }, inplace=True)

    # --- 4. Append results to the main CSV ---
    if results_csv_path.exists():
        print(f"Appending results to {results_csv_path.name}")
        existing_df = pd.read_csv(results_csv_path)
        final_df = pd.concat([existing_df, frame_df], ignore_index=True)
    else:
        print(f"Creating new results file: {results_csv_path.name}")
        final_df = frame_df

    final_df.to_csv(results_csv_path, index=False)
    print("Append complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Cellpose and analysis on a single image frame.")
    parser.add_argument("image_file", type=str, help="The full path to the image file to process.")
    args = parser.parse_args()

    image_to_process = Path(args.image_file)
    if not image_to_process.exists():
        raise FileNotFoundError(f"Input image not found: {image_to_process}")

    # Load the main project config
    # Assumes this script is in the 'scripts' directory
    config_path = Path(__file__).resolve().parents[1] / 'config.json'
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found at: {config_path}")

    with open(config_path, 'r') as f:
        project_config = json.load(f)

    process_frame(image_to_process, project_config)
