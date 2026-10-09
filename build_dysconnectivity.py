import argparse
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd


# ============================================================
# Load one pscalar file
# ============================================================

def load_pscalar(pscalar_file):
    """
    Load one parcellated CIFTI scalar file.

    Returns:
        values       - one disconnectivity value per parcel
        parcel_names - names of the parcels, if available
    """

    img = nib.load(pscalar_file)

    data = img.get_fdata()

    # A pscalar generally has shape:
    # (1, number_of_parcels)
    values = np.asarray(data).squeeze()

    if values.ndim != 1:
        raise ValueError(
            f"Expected one-dimensional parcel values in "
            f"{pscalar_file}, but got shape {data.shape}"
        )

    # The parcel axis is normally axis 1.
    parcel_axis = img.header.get_axis(1)

    parcel_names = list(parcel_axis.name)

    if len(parcel_names) != len(values):
        raise ValueError(
            f"Number of parcel names does not match number of "
            f"values in {pscalar_file}"
        )

    return values.astype(float), parcel_names


# ============================================================
# Extract subject ID from filename
# ============================================================

def get_subject_id(filename, prefix):
    """
    Extract subject ID from filenames such as:

        Glasser_1005653.pscalar.nii

    using prefix='Glasser_'.
    """

    name = filename.name

    expected_suffix = ".pscalar.nii"

    if not name.endswith(expected_suffix):
        raise ValueError(
            f"Unexpected filename: {name}"
        )

    stem = name[
        : -len(expected_suffix)
    ]

    if not stem.startswith(prefix):
        raise ValueError(
            f"Filename {name} does not begin with "
            f"expected prefix {prefix}"
        )

    subject_id = stem[
        len(prefix):
    ]

    return subject_id


# ============================================================
# Build subject × parcel matrix
# ============================================================

def build_matrix(
    input_dir,
    prefix,
):
    """
    Read all matching pscalar files and construct:

        rows    = subjects
        columns = Glasser parcels
        values  = disconnectivity
    """

    input_dir = Path(input_dir)

    pattern = (
        f"{prefix}*.pscalar.nii"
    )

    files = sorted(
        input_dir.glob(pattern)
    )

    if len(files) == 0:
        raise FileNotFoundError(
            f"No files matching {pattern} were found in "
            f"{input_dir}"
        )

    print(
        f"Found {len(files)} pscalar files"
    )

    rows = []

    reference_parcel_names = None

    for i, pscalar_file in enumerate(
        files,
        start=1,
    ):
        subject_id = get_subject_id(
            pscalar_file,
            prefix,
        )

        print(
            f"Loading {i}/{len(files)}: "
            f"{subject_id}"
        )

        values, parcel_names = load_pscalar(
            pscalar_file
        )

        # Use the first subject to define the expected parcel order.
        if reference_parcel_names is None:
            reference_parcel_names = parcel_names

            print(
                f"Number of parcels: "
                f"{len(reference_parcel_names)}"
            )

        else:
            if parcel_names != reference_parcel_names:
                raise ValueError(
                    f"Parcel ordering in {pscalar_file} does "
                    f"not match the first subject"
                )

        row = {
            "subject_id": subject_id
        }

        for parcel_name, value in zip(
            reference_parcel_names,
            values,
        ):
            row[parcel_name] = value

        rows.append(row)

    dataframe = pd.DataFrame(
        rows
    )

    dataframe = dataframe.sort_values(
        "subject_id"
    ).reset_index(
        drop=True
    )

    return dataframe


# ============================================================
# Quality-control summary
# ============================================================

def print_summary(dataframe):
    """
    Print simple checks before saving.
    """

    parcel_columns = [
        column
        for column in dataframe.columns
        if column != "subject_id"
    ]

    values = dataframe[
        parcel_columns
    ].to_numpy(
        dtype=float
    )

    print()
    print("==============================")
    print("Matrix summary")
    print("==============================")

    print(
        f"Subjects: {len(dataframe)}"
    )

    print(
        f"Parcels: {len(parcel_columns)}"
    )

    print(
        f"Matrix shape: {values.shape}"
    )

    print(
        f"Missing values: "
        f"{np.isnan(values).sum()}"
    )

    print(
        f"Minimum value: "
        f"{np.nanmin(values):.8g}"
    )

    print(
        f"Maximum value: "
        f"{np.nanmax(values):.8g}"
    )

    print(
        f"Mean value: "
        f"{np.nanmean(values):.8g}"
    )

    completely_zero = np.all(
        values == 0,
        axis=1,
    )

    print(
        f"Subjects with completely zero profile: "
        f"{completely_zero.sum()}"
    )


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Build a subject-by-parcel dysconnectivity matrix "
            "from Glasser pscalar files."
        )
    )

    parser.add_argument(
        "--input-dir",
        required=True,
        help=(
            "Directory containing Glasser pscalar files."
        ),
    )

    parser.add_argument(
        "--prefix",
        default="Glasser_",
        help=(
            "Filename prefix before the subject ID. "
            "Default: Glasser_"
        ),
    )

    parser.add_argument(
        "--output",
        required=True,
        help=(
            "Output CSV containing subjects as rows and "
            "parcels as columns."
        ),
    )

    args = parser.parse_args()

    dataframe = build_matrix(
        input_dir=args.input_dir,
        prefix=args.prefix,
    )

    print_summary(
        dataframe
    )

    output_path = Path(
        args.output
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    dataframe.to_csv(
        output_path,
        index=False,
    )

    print()
    print(
        f"Saved matrix to:"
    )
    print(
        output_path
    )


if __name__ == "__main__":
    main()