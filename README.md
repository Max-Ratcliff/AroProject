# Spore Germination Analysis Pipeline

This project contains the scripts and notebooks used to analyze spore germination from time-lapse microscopy images. The workflow consists of:

1.  **Segmentation:** `scripts/run_cellpose.py` - Segments spores from raw TIFF images.
2.  **Measurement:** `scripts/run_mask_analysis.py` - Measures spore properties using the generated masks.
3.  **Analysis:** `notebooks/MultiMeasure_Processing.ipynb` - Merges and analyzes data from different channels.

## How to Run

1.  Install dependencies: `pip install -r requirements.txt`
2.  Place raw data in the appropriate folders.
3.  Run the scripts in order.
