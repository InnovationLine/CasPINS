#!/usr/bin/env python
"""
Fix manuscript figures for Bioinformatics Advances submission.

Fixes:
1. Remove duplicate panel labels using precisely measured coordinates
2. Improve image quality using integer nearest-neighbor upscale (preserves
   screenshot pixel edges) followed by Lanczos adjustment to exact target size
3. Add clean, uniform panel labels (bold uppercase, white box, upper-left)
4. Save at 600 DPI TIFF with LZW compression

OUP Guidelines for combination figures: minimum 600 DPI, TIFF with LZW
"""

import os
import math
import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageFilter, ImageEnhance

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.join(SCRIPT_DIR, '..')
EXTRACTED_DIR = os.path.join(BENCHMARK_DIR, 'results', 'manuscript_figures', 'extracted')
OUT_DIR = os.path.join(BENCHMARK_DIR, 'results', 'manuscript_figures', 'publication')
PREVIEW_DIR = os.path.join(OUT_DIR, 'preview')
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(PREVIEW_DIR, exist_ok=True)

DPI = 600
PAGE_WIDTH_INCHES = 7.0
TARGET_WIDTH = int(PAGE_WIDTH_INCHES * DPI)

LABEL_FONT_SIZE_PX = 48

try:
    font_bold = ImageFont.truetype("arialbd.ttf", LABEL_FONT_SIZE_PX)
except Exception:
    try:
        font_bold = ImageFont.truetype("arial.ttf", LABEL_FONT_SIZE_PX)
    except Exception:
        font_bold = ImageFont.load_default()


# ---- Figure definitions with PRECISELY MEASURED label coordinates ----
# old_rect: (x1, y1, x2, y2) rectangle to white-out in ORIGINAL image coords
# new_pos:  (x, y) position for new clean label in ORIGINAL coords

FIGURES = {
    'Fig1_gRNA_Design': {
        'src': os.path.join(EXTRACTED_DIR, 'figure_6.jpeg'),
        # 1441x478, horizontal 2-panel, labels "A." top-left, "B." top-center-right
        'panels': [
            {'letter': 'A', 'old_rect': (0, 0, 52, 52), 'new_pos': (5, 3)},
            {'letter': 'B', 'old_rect': (695, 0, 750, 52), 'new_pos': (700, 3)},
        ],
    },
    'Fig2_TALEN_Design': {
        'src': os.path.join(EXTRACTED_DIR, 'figure_2.jpeg'),
        # 1056x1637, vertical 3-panel
        'panels': [
            {'letter': 'A', 'old_rect': (0, 0, 48, 52), 'new_pos': (5, 3)},
            {'letter': 'B', 'old_rect': (0, 530, 48, 585), 'new_pos': (5, 535)},
            {'letter': 'C', 'old_rect': (0, 1085, 52, 1145), 'new_pos': (5, 1090)},
        ],
    },
    'Fig3_Primer_Design_Input': {
        'src': os.path.join(EXTRACTED_DIR, 'figure_4.jpeg'),
        # 1442x1607, vertical 2-panel; B label confirmed at y=777-805 by pixel scan
        'panels': [
            {'letter': 'A', 'old_rect': (0, 0, 55, 60), 'new_pos': (5, 3)},
            {'letter': 'B', 'old_rect': (0, 770, 55, 815), 'new_pos': (5, 775)},
        ],
    },
    'Fig3_Primer_Design_Results': {
        'src': os.path.join(EXTRACTED_DIR, 'figure_7.jpeg'),
        # 1438x1603, vertical 2-panel (continuing C, D)
        'panels': [
            {'letter': 'C', 'old_rect': (0, 0, 48, 55), 'new_pos': (5, 3)},
            {'letter': 'D', 'old_rect': (0, 785, 50, 850), 'new_pos': (5, 790)},
        ],
    },
    'Fig4_Single_Indel': {
        'src': os.path.join(EXTRACTED_DIR, 'figure_5.jpeg'),
        # 1094x1599, vertical 3-panel
        'panels': [
            {'letter': 'A', 'old_rect': (0, 0, 48, 55), 'new_pos': (5, 3)},
            {'letter': 'B', 'old_rect': (0, 468, 45, 520), 'new_pos': (5, 473)},
            {'letter': 'C', 'old_rect': (0, 1015, 50, 1075), 'new_pos': (5, 1020)},
        ],
    },
    'Fig5_Batch_Analysis': {
        'src': os.path.join(EXTRACTED_DIR, 'figure_1.jpeg'),
        # 1206x1740, vertical 3-panel
        'panels': [
            {'letter': 'A', 'old_rect': (0, 0, 45, 50), 'new_pos': (5, 3)},
            {'letter': 'B', 'old_rect': (0, 615, 40, 665), 'new_pos': (5, 620)},
            {'letter': 'C', 'old_rect': (0, 1165, 55, 1228), 'new_pos': (5, 1170)},
        ],
    },
    'Graphical_Abstract': {
        'src': os.path.join(EXTRACTED_DIR, 'figure_3.jpg'),
        'panels': [],
    },
}


def sample_background_color(img_arr, rect, margin=5):
    """Sample the median background color from pixels adjacent to the rectangle."""
    x1, y1, x2, y2 = rect
    h, w = img_arr.shape[:2]

    samples = []
    # Right edge
    for y in range(max(0, y1), min(h, y2)):
        for dx in range(1, margin + 1):
            if x2 + dx < w:
                samples.append(img_arr[y, x2 + dx])
    # Bottom edge
    for x in range(max(0, x1), min(w, x2)):
        for dy in range(1, margin + 1):
            if y2 + dy < h:
                samples.append(img_arr[y2 + dy, x])
    # Top edge (above the rect)
    for x in range(max(0, x1), min(w, x2)):
        for dy in range(1, margin + 1):
            if y1 - dy >= 0:
                samples.append(img_arr[y1 - dy, x])

    if not samples:
        return (255, 255, 255)

    samples = np.array(samples)
    # Use the most common bright color (background is usually white/light gray)
    bright_mask = np.mean(samples, axis=1) > 180
    if np.sum(bright_mask) > 10:
        median_color = np.median(samples[bright_mask], axis=0).astype(int)
    else:
        median_color = np.median(samples, axis=0).astype(int)
    return tuple(median_color)


def smart_upscale(img, target_width):
    """
    Upscale screenshot using integer nearest-neighbor first, then Lanczos adjust.

    This preserves the crisp pixel edges of screenshots much better than
    direct Lanczos upscaling, which blurs everything.
    """
    orig_w, orig_h = img.size
    scale_needed = target_width / orig_w

    # Find the smallest integer scale >= scale_needed
    int_scale = math.ceil(scale_needed)
    # Cap at 4x to avoid excessive memory use
    int_scale = min(int_scale, 4)

    # Step 1: Integer nearest-neighbor scale (preserves pixel grid)
    int_w = orig_w * int_scale
    int_h = orig_h * int_scale
    img_nn = img.resize((int_w, int_h), Image.NEAREST)

    # Step 2: Lanczos to exact target (usually slight downscale which is high-quality)
    final_h = int(orig_h * (target_width / orig_w))
    img_final = img_nn.resize((target_width, final_h), Image.LANCZOS)

    return img_final


def process_figure(name, info):
    """Process one manuscript figure with quality fixes."""
    src = info['src']
    if not os.path.exists(src):
        print(f"  SKIP: {src} not found")
        return

    img = Image.open(src).convert('RGB')
    orig_w, orig_h = img.size
    img_arr = np.array(img)

    print(f"  Original: {orig_w}x{orig_h}")

    # Step 1: Remove old panel labels with precise white-out
    draw = ImageDraw.Draw(img)
    for panel in info['panels']:
        rect = panel['old_rect']
        bg_color = sample_background_color(img_arr, rect)
        # Expand rect by a few pixels to catch anti-aliasing edges
        x1 = max(0, rect[0] - 1)
        y1 = max(0, rect[1] - 1)
        x2 = min(orig_w, rect[2] + 3)
        y2 = min(orig_h, rect[3] + 3)
        draw.rectangle([x1, y1, x2, y2], fill=bg_color)

    # Step 2: Smart upscale using integer NN + Lanczos
    img = smart_upscale(img, TARGET_WIDTH)
    new_w, new_h = img.size
    scale = new_w / orig_w

    # Step 3: Apply sharpening for text clarity
    # Moderate unsharp mask to enhance edges without creating artifacts
    img = img.filter(ImageFilter.UnsharpMask(radius=1.2, percent=100, threshold=2))

    # Slight contrast boost
    enhancer = ImageEnhance.Contrast(img)
    img = enhancer.enhance(1.05)

    # Step 4: Add clean standardized panel labels
    draw = ImageDraw.Draw(img)

    for panel in info['panels']:
        letter = panel['letter']
        ox, oy = panel['new_pos']
        x = int(ox * scale)
        y = int(oy * scale)

        try:
            bbox = font_bold.getbbox(letter)
            tw = bbox[2] - bbox[0]
            th = bbox[3] - bbox[1]
            y_offset = bbox[1]
        except Exception:
            tw, th = LABEL_FONT_SIZE_PX, LABEL_FONT_SIZE_PX
            y_offset = 0

        pad = 5
        # White background box
        draw.rectangle(
            [x, y, x + tw + pad * 2, y + th + pad * 2],
            fill='white',
            outline=None
        )
        # Bold black letter
        draw.text((x + pad, y + pad - y_offset), letter, fill='black', font=font_bold)

    # Step 5: Save as 600 DPI TIFF with LZW compression
    out_path = os.path.join(OUT_DIR, f'{name}.tiff')
    img.save(out_path, format='TIFF', dpi=(DPI, DPI), compression='tiff_lzw')

    # Save PNG preview for quick visual check
    preview_w = 1800
    preview_scale = preview_w / new_w
    preview = img.resize((preview_w, int(new_h * preview_scale)), Image.LANCZOS)
    preview_path = os.path.join(PREVIEW_DIR, f'{name}.png')
    preview.save(preview_path, dpi=(150, 150))

    file_size_mb = os.path.getsize(out_path) / (1024 * 1024)
    print(f"  Output: {new_w}x{new_h} px @ {DPI} DPI ({file_size_mb:.1f} MB)")
    print(f"  TIFF:    {out_path}")
    print(f"  Preview: {preview_path}")


def main():
    print("=" * 70)
    print("MANUSCRIPT FIGURE FIX v2")
    print("- Remove duplicate panel labels (precise coordinates)")
    print("- Smart upscale: integer NN then Lanczos (sharp screenshots)")
    print("- Clean uniform panel labels: bold uppercase, white box")
    print("=" * 70)
    print(f"Settings: {DPI} DPI, {PAGE_WIDTH_INCHES}\" width = {TARGET_WIDTH}px")
    print()

    for name, info in FIGURES.items():
        print(f"Processing {name}...")
        process_figure(name, info)
        print()

    print("=" * 70)
    print("ALL FIGURES FIXED SUCCESSFULLY")
    print(f"Output directory: {OUT_DIR}")
    print(f"Preview directory: {PREVIEW_DIR}")
    print("=" * 70)


if __name__ == '__main__':
    main()
