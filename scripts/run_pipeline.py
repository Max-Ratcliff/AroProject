#!/usr/bin/env python3
"""
Spore Germination Analysis Pipeline
Orchestrates different processing modes for spore germination analysis.

Usage:
  - Batch mode (master mask):
    python run_pipeline.py batch
  - Single image mode (unique mask per frame):
    python run_pipeline.py single path/to/your/image.tif
  - Watch mode (for live acquisition):
    python run_pipeline.py watch --interval 300
"""

import sys
import subprocess
import time
from pathlib import Path
import json
from natsort import natsorted
import argparse


def run_script(script_path, description, args_list=None):
    """Executes a python script as a subprocess, handling arguments and errors."""
    start_time = time.time()
    print(f"\n🚀 Running: {description}...")
    command = [sys.executable, str(script_path)]
    if args_list:
        command.extend(args_list)

    try:
        subprocess.run(
            command,
            check=True,
            text=True,
            stdout=sys.stdout,
            stderr=subprocess.STDOUT,
        )
        duration = time.time() - start_time
        print(f"✅ Success ({duration:.1f}s)")
        return True
    except subprocess.CalledProcessError as e:
        duration = time.time() - start_time
        print(f"❌ FAILED: Command '{' '.join(e.cmd)}' returned non-zero exit status {e.returncode}. ({duration:.1f}s)")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Spore Germination Analysis Pipeline.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Available commands")

    # --- Batch Command ---
    parser_batch = subparsers.add_parser(
        "batch",
        help="Run the full batch processing pipeline (master mask) on a folder.\n"
             "Uses CURRENT_EXPERIMENT from config.json by default."
    )
    parser_batch.add_argument(
        "--input_dir", type=str,
        help="Path to the raw image folder. Overrides config.json."
    )
    parser_batch.add_argument(
        "--skip_cellpose", action="store_true",
        help="Skip Cellpose and use an existing master mask."
    )

    # --- Single Frame Command ---
    parser_single = subparsers.add_parser(
        "single",
        help="Process a single image file (generates unique mask and appends to CSV)."
    )
    parser_single.add_argument(
        "image_file", type=str,
        help="Path to the single .tif image to process."
    )

    # --- Watcher Command ---
    parser_watch = subparsers.add_parser(
        "watch",
        help="Watch a directory for new images and process them individually."
    )
    parser_watch.add_argument(
        "--interval", type=int, default=60,
        help="Time to wait between checking for new files, in seconds. Default: 60"
    )
    parser_watch.add_argument(
        "--watch_dir", type=str,
        help="Directory to watch. Overrides config.json."
    )

    args = parser.parse_args()
    scripts_dir = Path(__file__).parent

    # --- Execute Commands ---
    if args.command == "batch":
        # Load config to find all microscopy types for the experiment
        config_path = scripts_dir.parent / 'config.json'
        with open(config_path, 'r') as f:
            config = json.load(f)

        EXPERIMENT = config["CURRENT_EXPERIMENT"]
        try:
            microscopy_types = list(config["experiments"][EXPERIMENT].keys())
        except KeyError:
            print(f"Error: Experiment '{EXPERIMENT}' not found in config.json")
            return 1

        print(f"--- Starting Batch Mode for experiment '{EXPERIMENT}' ---")
        print(f"Found channels to process: {microscopy_types}")
        total_start_time = time.time()

        # Steps 1 & 2: Run Cellpose and Mask Analysis for each channel
        for m_type in microscopy_types:
            print(f"\n--- Processing channel: {m_type} ---")
            # Pass the microscopy type to the sub-scripts
            script_args = [f"--microscopy_type={m_type}"]

            if not args.skip_cellpose:
                if not run_script(scripts_dir / "run_cellpose.py", f"Cellpose ({m_type})", script_args):
                    return 1

            if not run_script(scripts_dir / "run_mask_analysis.py", f"Mask Analysis ({m_type})", script_args):
                return 1

        # Step 3: Run Multi-channel Processing
        # This script implicitly uses CURRENT_EXPERIMENT from config, which is correct.
        if len(microscopy_types) > 1:
            if not run_script(scripts_dir / "multimeasure_processing.py", "Multi-channel Processing"):
                return 1

        print(f"\n🎉 Batch pipeline completed successfully in {time.time() - total_start_time:.1f}s")

    elif args.command == "single":
        print("--- Starting Single Image Mode ---")
        if not run_script(
            scripts_dir / "process_single_frame.py",
            f"Processing {Path(args.image_file).name}",
            args_list=[args.image_file]
        ):
            return 1

    elif args.command == "watch":
        print(f"--- Starting Watch Mode (Interval: {args.interval}s) ---")
        print("Press Ctrl+C to stop.")

        # Determine which directory to watch
        if args.watch_dir:
            raw_data_dir = Path(args.watch_dir)
        else:
            # Load from config if not provided
            config_path = Path(__file__).resolve().parents[1] / 'config.json'
            with open(config_path, 'r') as f:
                config = json.load(f)
            DATA_ROOT = Path(config["DATA_ROOT"])
            EXPERIMENT = config["CURRENT_EXPERIMENT"]
            MICROSCOPY_TYPE = config["CURRENT_MICROSCOPY_TYPE"]
            exp_config = config["experiments"][EXPERIMENT][MICROSCOPY_TYPE]
            raw_data_dir = DATA_ROOT / "data" / "raw" / EXPERIMENT / exp_config["RAW_SUBFOLDER"]

        print(f"Watching for new .tif files in: {raw_data_dir}")
        processed_files = set()
        script_to_run = scripts_dir / "process_single_frame.py"

        while True:
            try:
                all_files = set(natsorted(raw_data_dir.glob("*.tif")))
                new_files = all_files - processed_files

                if new_files:
                    for file_path in natsorted(list(new_files)):
                        run_script(
                            script_to_run,
                            f"Processing new frame {file_path.name}",
                            args_list=[str(file_path)]
                        )
                        processed_files.add(file_path)
                else:
                    print("No new frames found. Waiting...", end="\r", flush=True)

                time.sleep(args.interval)
            except KeyboardInterrupt:
                print("\nWatch mode stopped by user.")
                break
            except FileNotFoundError:
                print(f"\nError: Watch directory not found: {raw_data_dir}")
                return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
