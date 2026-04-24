#!/usr/bin/env python3
"""
Build composite PNG images (all panels stitched into one figure) for
Supp Figs S22-S28 from the native-editable PPT-rendered multi-page PDFs.

Pipeline:
  1. For each `FigS{N}_*_native_editable.pdf`, run `pdftoppm` to rasterize
     each page to a high-res PNG (-r 300).
  2. Open per-page PNGs with PIL, crop whitespace, and stitch into a grid
     appropriate to the figure:
       S22 composite        : 7 pages (1 title + 6 panels) -> 2x3 grid (drop title)
       S23 pre/post/Δ heatmap: 1 page -> single PNG (just rename)
       S24 canonical sigs    : 7 pages (1 title + 6 panels) -> 2x3 grid
       S25 IBI vs IAE        : 3 pages (1 title + 2 panels) -> 1x2 grid
       S26 subject radar     : 4 pages (1 title + 3 panels) -> 1x3 grid
       S27 pre-spec primary  : 5 pages (1 title + 4 panels) -> 2x2 grid
       S28 platform concord. : 3 pages (1 title + 2 panels) -> 1x2 grid

  3. Add an A/B/C... panel letter overlay + a small super-title strip at
     the top of each composite, and save as PNG.

Output in manuscript/figures/*_composite.png alongside the existing
matplotlib composites (kept; these are the PPT-native composite variants).
"""
from __future__ import annotations
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path("/mnt/sda1/data/TNT/analysis/260424_nanostring")
FIG = ROOT / "manuscript" / "figures"

# Layout plans: {pdf_stem: (rows, cols, skip_leading_pages, supertitle)}
LAYOUTS = {
    "FigS22_NanoString_composite_native_editable":
        (2, 3, 1, "Supp Fig S22  NanoString PanCancer Immune orthogonal validation — 6-panel composite"),
    "FigS23_NanoString_prepostdelta_heatmap_native_editable":
        (1, 1, 0, "Supp Fig S23  NanoString pre / post / Δ composite heatmap"),
    "FigS24_NanoString_canonical_signatures_native_editable":
        (2, 3, 1, "Supp Fig S24  NanoString regulatory-grade signatures — pre vs post"),
    "FigS25_IBI_vs_IAE_fingerprint_native_editable":
        (1, 2, 1, "Supp Fig S25  Inflamed-but-Ineffective (IBI) vs Inflamed-Active-Effective (IAE)"),
    "FigS26_subject_radar_native_editable":
        (1, 3, 1, "Supp Fig S26  Subject-level deep-dive — s4, s2, s11"),
    "FigS27_NanoString_prespec_primary_paired_native_editable":
        (2, 2, 1, "Supp Fig S27  Pre-registered Arrow 5 rescue primary paired Δ"),
    "FigS28_NanoString_platform_concordance_native_editable":
        (1, 2, 1, "Supp Fig S28  NanoString vs RNA-seq platform concordance"),
}

# Panel letter per layout (read left-to-right, top-to-bottom)
PANEL_LETTERS = "ABCDEF"

DPI = 200  # rasterization DPI (300 too heavy for 30 slides; 200 is crisp for submission)


def rasterize_pdf(pdf_path: Path, out_dir: Path, dpi: int = DPI) -> list[Path]:
    """Run pdftoppm to convert all pages to PNG. Returns sorted list of PNG paths."""
    prefix = out_dir / pdf_path.stem
    subprocess.run(
        ["pdftoppm", "-png", "-r", str(dpi), str(pdf_path), str(prefix)],
        check=True, capture_output=True
    )
    return sorted(out_dir.glob(f"{pdf_path.stem}-*.png"))


def crop_whitespace(img: Image.Image, threshold: int = 252, pad: int = 20) -> Image.Image:
    """Crop outer white border (pixels >= threshold in all channels) and re-pad."""
    bbox_img = img.convert("L").point(lambda p: 0 if p < threshold else 255)
    bbox = bbox_img.getbbox()
    if bbox is None:
        return img
    # invert: we want regions where the image is NOT all-white
    inv = Image.eval(bbox_img, lambda p: 255 - p)
    bbox = inv.getbbox()
    if bbox is None:
        return img
    left, upper, right, lower = bbox
    left = max(0, left - pad); upper = max(0, upper - pad)
    right = min(img.width, right + pad); lower = min(img.height, lower + pad)
    return img.crop((left, upper, right, lower))


def overlay_letter(img: Image.Image, letter: str, font_size: int = 64) -> Image.Image:
    """Stamp a bold panel letter in the top-left corner."""
    out = img.copy()
    draw = ImageDraw.Draw(out)
    # Try a few font paths (Linux / DejaVuSans-Bold is common)
    font = None
    for candidate in ("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
                       "/usr/share/fonts/TTF/DejaVuSans-Bold.ttf"):
        if os.path.exists(candidate):
            try:
                font = ImageFont.truetype(candidate, font_size)
                break
            except Exception:
                pass
    if font is None:
        font = ImageFont.load_default()
    draw.text((24, 18), letter, fill=(17, 17, 17), font=font)
    return out


def stitch_grid(panel_images: list[Image.Image], rows: int, cols: int,
                 gap: int = 18, bg: tuple = (255, 255, 255)) -> Image.Image:
    """Stitch panel_images into a rows x cols grid; each panel padded to common size."""
    # normalise each panel to the same size (max width & max height across set)
    max_w = max(img.width for img in panel_images)
    max_h = max(img.height for img in panel_images)
    # center each on canvas of (max_w, max_h)
    normed = []
    for img in panel_images:
        canvas = Image.new("RGB", (max_w, max_h), bg)
        ox = (max_w - img.width) // 2
        oy = (max_h - img.height) // 2
        canvas.paste(img, (ox, oy))
        normed.append(canvas)
    grid_w = cols * max_w + (cols - 1) * gap
    grid_h = rows * max_h + (rows - 1) * gap
    out = Image.new("RGB", (grid_w, grid_h), bg)
    for idx, img in enumerate(normed):
        r = idx // cols
        c = idx % cols
        x = c * (max_w + gap)
        y = r * (max_h + gap)
        out.paste(img, (x, y))
    return out


def add_supertitle(img: Image.Image, title: str, font_size: int = 42,
                    pad: int = 30, bg: tuple = (255, 255, 255)) -> Image.Image:
    """Prepend a white strip with the supertitle."""
    font = None
    for candidate in ("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
                       "/usr/share/fonts/TTF/DejaVuSans-Bold.ttf"):
        if os.path.exists(candidate):
            try:
                font = ImageFont.truetype(candidate, font_size)
                break
            except Exception:
                pass
    if font is None:
        font = ImageFont.load_default()
    strip_h = font_size + 2 * pad
    out = Image.new("RGB", (img.width, img.height + strip_h), bg)
    draw = ImageDraw.Draw(out)
    # measure
    try:
        tb = draw.textbbox((0, 0), title, font=font)
        tw, th = tb[2] - tb[0], tb[3] - tb[1]
    except Exception:
        tw, th = len(title) * (font_size // 2), font_size
    draw.text(((img.width - tw) // 2, pad), title, fill=(17, 17, 17), font=font)
    out.paste(img, (0, strip_h))
    return out


def build_composite(pdf_stem: str, rows: int, cols: int,
                     skip_leading: int, title: str, tmpdir: Path) -> Path:
    pdf_path = FIG / f"{pdf_stem}.pdf"
    if not pdf_path.exists():
        print(f"SKIP (no pdf): {pdf_path}")
        return None
    png_paths = rasterize_pdf(pdf_path, tmpdir)
    # drop header slide(s) if requested
    panel_pngs = png_paths[skip_leading:]
    expected = rows * cols
    if len(panel_pngs) != expected:
        print(f"WARN {pdf_stem}: expected {expected} panels, got {len(panel_pngs)} (rows={rows}, cols={cols})")
        # if more than expected, truncate; if fewer, pad with blanks
        if len(panel_pngs) > expected:
            panel_pngs = panel_pngs[:expected]
    # Load, crop whitespace, add panel letter
    imgs = []
    for i, p in enumerate(panel_pngs):
        img = Image.open(p).convert("RGB")
        img = crop_whitespace(img, pad=25)
        if len(panel_pngs) > 1:  # single-panel figs get no letter
            img = overlay_letter(img, PANEL_LETTERS[i])
        imgs.append(img)
    # pad to grid length if short
    while len(imgs) < expected:
        canvas = Image.new("RGB", imgs[0].size, (255, 255, 255))
        imgs.append(canvas)
    grid = stitch_grid(imgs, rows, cols)
    with_title = add_supertitle(grid, title)
    out_path = FIG / f"{pdf_stem.replace('_native_editable', '')}_composite_stitched.png"
    # resize so longer side <= 4000 px for a manageable file
    max_dim = 4000
    if max(with_title.size) > max_dim:
        scale = max_dim / max(with_title.size)
        new_size = (int(with_title.width * scale), int(with_title.height * scale))
        with_title = with_title.resize(new_size, Image.LANCZOS)
    with_title.save(out_path, "PNG", optimize=True)
    print(f"wrote {out_path.name}  ({with_title.width}x{with_title.height} px, "
          f"{out_path.stat().st_size // 1024} KB)")
    return out_path


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        tmp_p = Path(tmp)
        for stem, (rows, cols, skip, title) in LAYOUTS.items():
            build_composite(stem, rows, cols, skip, title, tmp_p)
    print("\n=== All composite PNGs written to", FIG, "===")
