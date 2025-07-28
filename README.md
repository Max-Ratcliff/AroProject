# Spore Germination Analysis Pipeline

This project contains the scripts and notebooks used to analyze spore germination from time-lapse microscopy images. The workflow consists of:

1. **Setup** - modify .env file to set your root directory this is my directory structure:
    .
    ├── data
    │   ├── processed
    │   │   ├── cellpose
    │   │   │   ├── Leticia_M4576_s2_PH
    │   │   │   ├── Leticia_M4576_s2_ThT
    │   │   │   └── Pulses_ThT_ThT
    │   │   ├── Leticia_M4576_s2_Fiji
    │   │   │   └── Processed_Data
    │   │   └── Pulses_ThT_Fiji
    │   │       └── Processed_Data
    │   └── raw
    │       ├── autofluorescence
    │       ├── Leticia_M4576_s2
    │       │   ├── M4576_s2_PH
    │       │   │   └── M4576_s2_PH_masks
    │       │   └── M4576_s2_ThT
    │       └── Pulses_ThT
    │           ├── M4581_s2_PH_stabilized
    │           ├── M4581_s2_ThT_stabilized
    │           ├── M4584_s1_PH_stabilized
    │           ├── M4584_s1_ThT_stabilized
    │           ├── M4867_s2_PH_stabilized
    │           ├── M4867_s2_ThT_stabilized
    │           ├── M4867_s3_PH_stabilized
    │           ├── M4867_s3_ThT_stabilized
    │           ├── M4868_s5_PH_stabilized
    │           ├── M4868_s5_ThT_stabilized
    │           ├── M4868_s6_PH_stabilized
    │           ├── M4868_s6_ThT_stabilized
    │           ├── M6881_s2_PH_stabilized
    │           ├── M6881_s2_ThT_stabilized
    │           ├── M6881_s5_PH_stabilized
    │           ├── M6881_s5_ThT_stabilized
    │           ├── M6881_s6_PH_stabilized
    │           ├── M6881_s6_ThT_stabilized
    │           └── M6881_s6_w1PH80_stabilized
    ├── notebooks
    ├── reports
    │   └── figures
    └── scripts
2. **Segmentation:** `scripts/run_cellpose.py` - Segments spores from raw TIFF images.
3. **Measurement:** `scripts/run_mask_analysis.py` - Measures spore properties using the generated masks.
4. **Analysis:** `notebooks/MultiMeasure_Processing.ipynb` - Merges and analyzes data from different channels.
5. **Metrics Calculation:** `notebooks/exp_metrics.ipynb` - Calculates experimental metrics based on the processed data.

## How to Run

1. Install dependencies: `pip install -r requirements.txt`
2. Place raw data in the appropriate folders.
3. Run the scripts in order.
