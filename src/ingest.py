# list_series_sitk.py
from pathlib import Path
import SimpleITK as sitk
import pydicom

def list_series(dicom_dir: str):
    dicom_dir = Path(dicom_dir)
    r = sitk.ImageSeriesReader()
    sids = r.GetGDCMSeriesIDs(str(dicom_dir))
    if not sids:
        raise SystemExit("No DICOM series found")

    rows = []
    for sid in sids:
        files = r.GetGDCMSeriesFileNames(str(dicom_dir), sid)
        rows.append((sid, len(files)))

    # sort by slice count desc
    rows.sort(key=lambda x: x[1], reverse=True)

    print("Found series:")
    for sid, n in rows:
        print(f"- SeriesInstanceUID={sid} | files={n}")

    best_sid, best_n = rows[0]
    print("\nSelected (largest):")
    print(f"SeriesInstanceUID={best_sid} | files={best_n}")
    # If you want the file list to pass to the next step:
    best_files = r.GetGDCMSeriesFileNames(str(dicom_dir), best_sid)
    return best_sid, best_files

def sorted_files_sitk(dicom_dir: str):
    r = sitk.ImageSeriesReader()
    sids = r.GetGDCMSeriesIDs(dicom_dir)
    if not sids: raise RuntimeError("No series found")
    sid = sids[0]  # or choose the largest series first
    files = r.GetGDCMSeriesFileNames(dicom_dir, sid)  # <- sorted
    return sid, files

def convert_hounsfield_unit(dicom_dir: str):
    i, files = sorted_files_sitk(dicom_dir)
    for f in files:
        dataset = pydicom.dcmread("your_dicom_file.dcm")
        raw_pixel_data = dataset.pixel_array
        # HU represents tissue density or radiodensity
        hounsfield_units = pydicom.pixel_data_handlers.util.apply_modality_lut(raw_pixel_data, dataset)
    return hounsfield_units

if __name__ == "__main__":
    # change to your folder path
    list_series("../data/CASE1")