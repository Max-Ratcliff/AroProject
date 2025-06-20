from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from cellpose import core, io, models  # , plot
from natsort import natsorted
from dotenv import load_dotenv
import os
# from tqdm import trange

io.logger_setup()

if core.use_gpu() is False:
    raise ImportError("no GPU access")

model = models.CellposeModel(gpu=True)

load_dotenv()
ROOT = os.getenv("DATA_ROOT")
if ROOT is None:
    raise ValueError("DATA_ROOT environment variable is not set. Please set it to the root directory of your data.")

# EXPERIMENT = "Leticia_M4576_s2"
EXPERIMENT = "Pulses_ThT"
CHANNEL = "ThT"
ROOT = Path(ROOT)
# DATA_DIR = ROOT / f"data/raw/{EXPERIMENT}/M4576_s2_PH/"
DATA_DIR = ROOT / f"data/raw/{EXPERIMENT}/M4581_s2_PH_stabilized/"

OUTPUT_DIR = ROOT / f"data/processed/cellpose/{EXPERIMENT}_{CHANNEL}"

OUTPUT_DIR = Path(OUTPUT_DIR)
if not OUTPUT_DIR.exists():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Created output directory: {OUTPUT_DIR}")

path = Path(DATA_DIR)

image_ext = ".tif"

files = natsorted(
    [
        file
        for file in path.glob("*" + image_ext)
        if "_masks" not in file.name and "_flows" not in file.name
    ]
)

if len(files) == 0:
    raise FileNotFoundError(f"No files found with extension {image_ext} in {path}")
else:
    print(f"Found {len(files)} files with extension {image_ext} in {path}")

for file in files:
    print(f"{file.name}")

img = io.imread(files[0])  # read first image to get shape
print(f"Image shape: {img.shape}")
# if i remember right we need to add a channel dimension because our tifs onlyh have one channel
# if len(img.shape) == 2:
#     img = img[:, :, np.newaxis]

# in the note book they do some weird channel logic ill paste it bellow in a comment i dont think we need it though
# selected_channels = []
# for i, c in enumerate([first_channel, second_channel, third_channel]):
#   if c == 'None':
#     continue
#   if int(c) > img.shape[-1]:
#     assert False, 'invalid channel index, must have index greater or equal to the number of channels'
#   if c != 'None':
#     selected_channels.append(int(c))
#

img_selected_channels = np.zeros_like(img)
img_selected_channels[:, :] = img[:, :]

imgs = [img]

flow_threshold = 0.4
cellprob_threshold = 0.0
tile_norm_blocksize = 0

fig = plt.figure(figsize=(12, 5))
plt.imshow(img)
# plt.title('Original Image')
plt.axis('off')
plt.tight_layout()
# plt.show()
fig.savefig(OUTPUT_DIR / "original_image.jpg", bbox_inches='tight', pad_inches=0)
print(f"Original image saved to {OUTPUT_DIR / 'original_image.jpg'}")

print("Running Cellpose model...")
print(f"Image shape: {img.shape}")
print(f"Flow threshold: {flow_threshold}")
print(f"Cell probability threshold: {cellprob_threshold}")

masks, flows, styles = model.eval(
    imgs,
    flow_threshold=flow_threshold,
    cellprob_threshold=cellprob_threshold,
)

print("Model evaluation complete.")
print(f"Number of masks found: {len(masks)}")
print(f"Number of flows found: {len(flows)}")
print(f"Number of styles found: {len(styles)}")

print(f"Mask shape: {masks[0].shape}")
# print(f"Flow shape: {flows.shape}")
# Save masks and flows
output_file_path = OUTPUT_DIR / files[0].name

# The save_masks function expects lists of images, masks, flows, and filenames
io.save_masks(
    [img],               # Original image, in a list
    [masks[0]],          # Masks, in a list
    [flows[0]],          # Flows, in a list
    [str(output_file_path)]  # File path, in a list
)
print(f"Masks and flows saved to {OUTPUT_DIR}")

# fig = plt.figure(figsize=(12, 5))
# plot.show_segmentation(fig, img, masks[0], flows[0])
# plt.tight_layout()
# plt.show()
# fig.savefig(test_out / "segmentation_result.jpg", bbox_inches='tight', pad_inches=0)

# Assuming 'flows' is the output from model.eval() for a single image
flow_components = flows[0]

# Unpack the components
rgb_flow_image = flow_components[0]
dP = flow_components[1]
cellprob = flow_components[2]

# Create subplots to display them
fig, axs = plt.subplots(1, 3, figsize=(18, 6))

# 1. Cell Probability
axs[0].imshow(cellprob, cmap='viridis')
axs[0].set_title('Cell Probability')
axs[0].axis('off')

# 2. Raw dP[0] (Y-flow) - less intuitive
axs[1].imshow(dP[0], cmap='coolwarm')
axs[1].set_title('Raw Flow (Y-component)')
axs[1].axis('off')

# 3. RGB Flow Image - Most useful view
axs[2].imshow(rgb_flow_image)
axs[2].set_title('RGB Flow Visualization')
axs[2].axis('off')

plt.tight_layout()

fig.savefig(OUTPUT_DIR / "flow_visualization.jpg", bbox_inches='tight', pad_inches=0)
print(f"Flow visualization saved to {OUTPUT_DIR / 'flow_visualization.jpg'}")
