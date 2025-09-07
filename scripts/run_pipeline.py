#!/usr/bin/env python3
"""
Sppoer Germination Analysis Pipeline
Runs all analysis scripts in order with timing and clear output.

Usage: python run_pipeline.py [step1|step2|...]
"""

import sys
import subprocess
import time
from pathlib import Path


def run_script(script_path, description):
    """Execute script and return (success, duration)"""
    start = time.time()
    print(f"🚀 {description}")
    result = subprocess.run([sys.executable, str(script_path)], stdout=sys.stdout, stderr=sys.stdout, text=True)
    return result.returncode == 0, time.time() - start


def main():
    scripts = [
        ("run_cellpose.py", "Cellpose Segmentation"),
        ("run_mask_analysis.py", "Mask Analysis"),
        ("multimeasure_processing.py", "Multi-channel Processing"),
        ("exp_metrics.py", "Experiment Metrics")
    ]
    start = time.time()
    for script, desc in scripts:
        success, duration = run_script(Path(__file__).parent / script, desc)
        print(f"✅ {desc} - {duration:.1f}s")
    print(f"✅ Total - {time.time() - start:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
