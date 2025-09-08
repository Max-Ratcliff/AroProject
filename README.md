# Spore Germination Analysis Pipeline

This project contains a suite of Python scripts to analyze spore germination from time-lapse microscopy images. The pipeline is designed to be flexible, supporting both post-experiment batch processing and live analysis as images are acquired.

## Directory Structure

The project expects the following directory structure:

```
.
├── config.json
├── data
│   ├── processed
│   │   ├── cellpose
│   │   │   └── Leticia_M4576_s2_PhC
│   │   ├── Leticia_M4576_s2_Fiji
│   │   │   └── Processed_Data
│   │   └── ...
│   └── raw
│       ├── Leticia_M4576_s2
│       │   ├── M4576_s2_PH
│       │   └── M4576_s2_ThT
│       └── ...
├── notebooks
├── reports
│   └── figures
├── requirements.txt
└── scripts
    ├── run_pipeline.py
    ├── run_cellpose.py
    ├── run_mask_analysis.py
    ├── multimeasure_processing.py
    ├── exp_metrics.py
    └── process_single_frame.py
```

## How to Run

### 1. Installation

First, ensure you have a compatible Python environment (e.g., using Conda). Then, install the required dependencies:

```bash
pip install -r requirements.txt
```

### 2. Configuration

All pipeline settings are managed in the `config.json` file. Before running, you must:

1.  **Set `DATA_ROOT`**: Update this to the absolute path of the `AroProject` directory on your machine.
2.  **Define Experiments**: Add your experiments to the `experiments` object. For each experiment, define the microscopy channels (e.g., `PhC`, `ThT`) and specify their `RAW_SUBFOLDER` and the desired `MASK_FILENAME`.
3.  **Set Current Task**: To run the pipeline, set `CURRENT_EXPERIMENT` and `CURRENT_MICROSCOPY_TYPE` to point to the data you want to process.

### 3. Running the Pipeline

The main entry point is the `run_pipeline.py` script, which offers several modes of operation.

#### Batch Mode

This mode processes an entire experiment at once. It automatically finds all microscopy channels (e.g., `PhC` and `ThT`) for the `CURRENT_EXPERIMENT` set in `config.json`, runs segmentation and analysis on each, and then merges the results for final processing and plotting.

```bash
# Run the full pipeline on the experiment defined in config.json
python scripts/run_pipeline.py batch
```

#### Single Image Mode

This mode is for processing a single image file. It generates a unique mask for that frame, measures its properties, and appends the results to the corresponding experiment's CSV file. This is ideal for live analysis where images arrive one at a time.

```bash
# Process a single image file
python scripts/run_pipeline.py single /path/to/your/image.tif
```

#### Watch Mode

This mode continuously monitors a directory for new images. When a new `.tif` file appears, it automatically triggers the "single image" processing for that file. This is perfect for automating analysis during live microscope acquisition.

```bash
# Watch the configured experiment directory, checking every 5 minutes (300s)
python scripts/run_pipeline.py watch --interval 300
```
