from pathlib import Path
import numpy as np
import pydicom
from pydicom.pixel_data_handlers.util import apply_modality_lut
import argparse, json, re, hashlib
import json
import matplotlib.pyplot as plt

def retrieve_dataset(dataset_dir: str):
    datasets = []
    data = Path(dataset_dir)

    count = 0
    for folder in data.iterdir():
        if folder.is_dir():
            datasets.append(folder)
            print(f"Dataset #{count}: {folder.name}")
            count += 1
    choice = input("Which dataset would you like to process? (enter number) ")

    try:
        idx = int(choice)
        return datasets[idx]
    except (ValueError, IndexError):
        print("Invalid choice.")
        return None

def dicom_directory(dataset_dir: str):
    """
    Returns { study_uid: { series_uid: { "desc": str, "modality": str, "files": [paths...] } } }
    Recurses through the dataset folder; you do NOT need to know its internal tree.
    """
    index: dict[str, dict[str, dict]] = {}  # {study_uid: {series_uid: {...}}}
    root = Path(dataset_dir)
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        try:
            ds = pydicom.dcmread(str(path), stop_before_pixels=True)
        except Exception:
            continue

        study_uid  = getattr(ds, "StudyInstanceUID",  None)
        series_uid = getattr(ds, "SeriesInstanceUID", None)
        if not (study_uid and series_uid):
            continue

        series_map = index.setdefault(study_uid, {})
        if series_uid not in series_map:
            series_map[series_uid] = {
                "desc": getattr(ds, "SeriesDescription", "Series"),
                "modality": getattr(ds, "Modality", "NA"),
                "files": []
            }
        series_map[series_uid]["files"].append(str(path))
    return index


def convert_files_to_hu(files: list[str]):
    # sort by InstanceNumber, then ImagePositionPatient[2], then name
    def sort_key(path):
        ds = pydicom.dcmread(path, stop_before_pixels=True)
        inst = getattr(ds, "InstanceNumber", None)
        if inst is not None: return (0, int(inst))
        ipp = getattr(ds, "ImagePositionPatient", None)
        if ipp is not None:  return (1, float(ipp[2]))
        return (2, Path(path).name)

    files_sorted = sorted(files, key=sort_key)

    slices, ipps = [], []
    ds0 = None
    for index, f in enumerate(files_sorted):
        ds = pydicom.dcmread(f)
        if index == 0: ds0 = ds
        hu = apply_modality_lut(ds.pixel_array, ds).astype(np.float32)
        slices.append(hu)
        ipp = getattr(ds, "ImagePositionPatient", [0, 0, index])
        ipps.append(np.array(ipp, dtype=float))

    vol_hu = np.stack(slices, axis=0)  # (Z,Y,X)
    row_mm, col_mm = map(float, ds0.PixelSpacing)
    dz = np.diff([p[2] for p in ipps])
    z_mm = float(np.median(np.abs(dz))) if len(dz) else float(getattr(ds0, "SliceThickness", 1.0))
    spacing_zyx = (z_mm, row_mm, col_mm)

    # Anatomy of a DICOM file!
    meta = {
        "StudyInstanceUID": getattr(ds0, "StudyInstanceUID",  "NA"),
        "SeriesInstanceUID": getattr(ds0, "SeriesInstanceUID", "NA"),
        "SeriesDescription": getattr(ds0, "SeriesDescription", "Series"),
        "Modality": getattr(ds0, "Modality", "NA"),
        "PhotometricInterpretation": getattr(ds0, "PhotometricInterpretation", "MONOCHROME2"),
        "NumSlices": int(vol_hu.shape[0]),
        "PixelSpacing": [row_mm, col_mm],
        "SliceSpacingEstimated": z_mm,
        "Window": {"center": 50, "width": 350}, # Clarification: this changes the window for the METADATA, not the actual volume.npy file
    }
    return vol_hu, spacing_zyx, meta


def window_to_unit(arr_hu: np.ndarray, center=50.0, width=350.0, invert_for_display=False):
    #6: Normalize to [0,1] (for volumes)
    lo, hi = center - width/2, center + width/2
    arr = np.clip(arr_hu, lo, hi)
    arr01 = (arr - lo) / (hi - lo + 1e-6)
    if invert_for_display:
        arr01 = 1.0 - arr01
    return arr01.astype(np.float32)


# Main Function: Get data, push it out
def main():
    # # Gets dataset (data/[dataset_name]) via a command line argument
    # ap = argparse.ArgumentParser()
    # ap.add_argument("--dataset", required=True, help="Path to one dataset folder under data/ (e.g., data/RSNA2022)")
    # args = ap.parse_args()

    dataset_dir = retrieve_dataset("E:\Spine-Mets-CT-data") # Dataset folder is on a freaking usb .-.
    if dataset_dir is None:
        print("No dataset selected.")
        return
    else:
        print("Chosen dataset:", dataset_dir)

    dataset_name = dataset_dir.name

    idx = dicom_directory(dataset_dir)
    if not idx:
        raise SystemExit(f"No studies found under {dataset_dir}")

    for study_uid, series_map in idx.items():
        #1: Largest Series in the Study
        series_uid, entry = max(series_map.items(), key=lambda kv: len(kv[1]["files"]))
        files = entry["files"]

        # 2-5 (the bulk of the ingestion, more so extraction, occurs here)
        vol_hu, spacing, meta = convert_files_to_hu(files)
        #6: Windowing and Normalizing vol_hu (the original vol_hu is now not needed)
        vol01 = window_to_unit(vol_hu, center=50, width=350)

        # Moving Out!
        #I: Assertions
        # Checking if volume.npy was actually written into
        vol = np.load(case_dir / "volume.npy", mmap_mode="r")
        assert vol.dtype == np.float32 and vol.ndim == 3 and vol.shape[0] >= 10
        # Checking if spacing is correct and all three decimals/floats are greater than 0
        s = json.loads((case_dir / "spacing.json").read_text())["spacing_zyx"]
        assert len(s) == 3 and all(float(x) > 0 for x in s)
        # Checking for failed save paths
        assert (case_dir / "preview.png").stat().st_size > 0
        
        #II: Getting Series, Study, and Series Desc for filenames and organization
        series_uid  = meta.get("SeriesInstanceUID", "NA")
        study_uid   = meta.get("StudyInstanceUID",  "NA")      # ensure convert_to_hu records this
        series_desc = meta.get("SeriesDescription", "Series")  # ensure convert_to_hu records this

        #III: Getting Preview Image
        z = vol01.shape[0] // 2
        preview_img = vol01[z].copy()
        if meta.get("PhotometricInterpretation","MONOCHROME2") == "MONOCHROME1":
            preview_img = 1.0 - preview_img
        # plt.imsave(out / "preview.png", preview_img, cmap="gray")

        #IV: Hashing UIDs for my sanity (UIDs are loooong)
        def short_id(uid: str, n=8): return hashlib.sha1(uid.encode()).hexdigest()[:n]
        def safe(s: str, n=24): return re.sub(r'[^A-Za-z0-9._-]+', '_', s or "NA")[:n]

        # Folder name: SER_<hash(series)>__STU_<hash(study)>__<sanitized series description>
        case_dir = Path("out") / dataset_name / f"STDY_{short_id(study_uid)}__SRS_{short_id(series_uid)}__{safe(series_desc)}" / "ingest"
        case_dir.mkdir(parents=True, exist_ok=True)

        save_ingestion(case_dir, vol01, spacing, meta, preview_img)
        

def save_ingestion(out_dir, vol01, spacing_zyx, meta, preview_img=None):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Save Series, Study, Metadata, and Preview Image
    np.save(out / "volume.npy", vol01.astype(np.float32))
    (out / "spacing.json").write_text(json.dumps({"spacing_zyx": spacing_zyx}, indent=2))
    (out / "meta.json").write_text(json.dumps(meta, indent=2))
    if preview_img.ndim != 2:
        return
    plt.imsave(out / "preview.png", preview_img, cmap="gray")