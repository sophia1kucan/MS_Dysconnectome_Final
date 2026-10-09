import argparse
import gzip
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, save_npz


# Load files

def load_dot(path):

    return np.loadtxt(path)


def load_dot_gz(path):

    with gzip.open(path, "rt") as file:
        return np.loadtxt(file)


def load_waytotal(path):

    with open(path, "r") as file:
        value = float(file.read().strip())

    if value <= 0:
        raise ValueError(
            f"Waytotal must be greater than zero: {path}"
        )

    return value


def load_subject(subject_dir, matrix_folder):
    """
    Load coordinates, connectivity matrices and waytotal values
    for one subject.

    Expected structure:

    subject_dir/
        mat2/
            omatrix2_L.white.resampled/
            omatrix2_R.white.resampled/

        mat2_pd/
            omatrix2_L.white.resampled/
            omatrix2_R.white.resampled/
    """

    subject_dir = Path(subject_dir)

    left_folder = (
        subject_dir
        / matrix_folder
        / "omatrix2_L.white.resampled"
    )

    right_folder = (
        subject_dir
        / matrix_folder
        / "omatrix2_R.white.resampled"
    )

    left_coords_path = (
        left_folder
        / "tract_space_coords_for_fdt_matrix2"
    )

    right_coords_path = (
        right_folder
        / "tract_space_coords_for_fdt_matrix2"
    )

    left_connectivity_path = (
        left_folder
        / "fdt_matrix2.dot.gz"
    )

    right_connectivity_path = (
        right_folder
        / "fdt_matrix2.dot.gz"
    )

    left_waytotal_path = (
        left_folder
        / "waytotal"
    )

    right_waytotal_path = (
        right_folder
        / "waytotal"
    )

    required_files = [
        left_coords_path,
        right_coords_path,
        left_connectivity_path,
        right_connectivity_path,
        left_waytotal_path,
        right_waytotal_path,
    ]

    missing_files = [
        path
        for path in required_files
        if not path.exists()
    ]

    if missing_files:
        missing_text = "\n".join(
            str(path)
            for path in missing_files
        )

        raise FileNotFoundError(
            f"Subject {subject_dir.name} is missing files:\n"
            f"{missing_text}"
        )

    coords_L = load_dot(left_coords_path)
    coords_R = load_dot(right_coords_path)

    conn_L = load_dot_gz(left_connectivity_path)
    conn_R = load_dot_gz(right_connectivity_path)

    waytotal_L = load_waytotal(left_waytotal_path)
    waytotal_R = load_waytotal(right_waytotal_path)

    return (
        coords_L,
        coords_R,
        conn_L,
        conn_R,
        waytotal_L,
        waytotal_R,
    )


# Coordinate handling

def coord_key(coord):
    #Convert one coordinate row into a rounded tuple.


    return tuple(
        np.round(
            np.asarray(coord, dtype=float)
        )
    )


def build_union_coordinate_space(all_coordinates):
    #Build a shared coordinate space containing every unique
    #tract-space coordinate found across all subjects.
  

    union = set()

    for coordinates in all_coordinates:
        for coordinate in coordinates:
            union.add(coord_key(coordinate))

    union_rows = sorted(union)

    coordinate_to_index = {
        coordinate: index
        for index, coordinate in enumerate(union_rows)
    }

    union_coordinates = np.asarray(
        union_rows,
        dtype=float
    )

    return (
        union_coordinates,
        coordinate_to_index,
    )


def build_subject_mappings(
    all_coordinates,
    coordinate_to_index,
):
    #Build one mapping per subject from the subject's local
    #tract-space coordinate index into the shared union space.

    subject_mappings = []

    for coordinates in all_coordinates:
        mapping = np.asarray(
            [
                coordinate_to_index[
                    coord_key(coordinate)
                ]
                for coordinate in coordinates
            ],
            dtype=np.int64,
        )

        subject_mappings.append(mapping)

    return subject_mappings


# Connectivity 

def project_subject(
    connectivity,
    subject_to_union,
    n_union,
):
    """
    Project one subject's sparse connectivity matrix into the
    shared union-coordinate space.

    The connectivity file is expected to contain three columns:

        seed_index  target_index  weight

    FSL indices are one-based, so they are converted to zero-based
    Python indices.
    """

    connectivity = np.asarray(connectivity)

    if connectivity.ndim != 2:
        raise ValueError(
            "Connectivity data must be a two-dimensional array"
        )

    if connectivity.shape[1] < 3:
        raise ValueError(
            "Connectivity data must contain at least three columns"
        )

    if connectivity.shape[0] == 0:
        raise ValueError(
            "Connectivity data contains no entries"
        )

    seed = (
        connectivity[:, 0].astype(np.int64)
        - 1
    )

    target = (
        connectivity[:, 1].astype(np.int64)
        - 1
    )

    weights = connectivity[:, 2].astype(
        np.float64
    )

    if seed.min() < 0:
        raise ValueError(
            "Seed indices must be positive before conversion "
            "from one-based indexing"
        )

    if target.min() < 0:
        raise ValueError(
            "Target indices must be positive before conversion "
            "from one-based indexing"
        )

    if target.max() >= len(subject_to_union):
        raise ValueError(
            "Target indices exceed the subject tract-space "
            "coordinate table"
        )

    target_union = subject_to_union[target]

    n_seed = int(seed.max()) + 1

    projected = coo_matrix(
        (
            weights,
            (
                seed,
                target_union,
            ),
        ),
        shape=(
            n_seed,
            n_union,
        ),
    ).tocsr()

    return projected


def apply_waytotal_normalisation(
    matrix,
    waytotal,
):
    """
    Divide every connectivity weight by the subject's waytotal.
    """

    if waytotal <= 0:
        raise ValueError(
            "Waytotal must be greater than zero"
        )

    matrix = matrix.astype(
        np.float64
    )

    return matrix / float(waytotal)


def row_normalise(matrix):
    """
    Normalise every cortical seed row so that its non-zero
    connectivity values sum to one.
    """

    row_sums = np.asarray(
        matrix.sum(axis=1)
    ).ravel()

    inverse_row_sums = np.zeros(
        row_sums.shape,
        dtype=np.float64,
    )

    nonzero_rows = row_sums != 0

    inverse_row_sums[nonzero_rows] = (
        1.0
        / row_sums[nonzero_rows]
    )

    return matrix.multiply(
        inverse_row_sums[:, np.newaxis]
    ).tocsr()


def build_subject_matrices(
    connectivity,
    subject_to_union,
    n_union,
    waytotal,
):
    #Need both normalised matrix versions for one subject.


    projected = project_subject(
        connectivity=connectivity,
        subject_to_union=subject_to_union,
        n_union=n_union,
    )

    waytotal_matrix = (
        apply_waytotal_normalisation(
            matrix=projected,
            waytotal=waytotal,
        )
    )

    row_matrix = row_normalise(
        waytotal_matrix
    )

    return (
        waytotal_matrix,
        row_matrix,
    )


# Hemisphere processing

def process_hemisphere(
    all_coordinates,
    all_connectivity,
    all_waytotals,
):


    n_subjects = len(all_coordinates)

    if n_subjects == 0:
        raise ValueError(
            "No subjects were supplied"
        )

    if len(all_connectivity) != n_subjects:
        raise ValueError(
            "Number of connectivity matrices does not match "
            "number of coordinate arrays"
        )

    if len(all_waytotals) != n_subjects:
        raise ValueError(
            "Number of waytotal values does not match "
            "number of subjects"
        )

    print("Building union coordinate space...")

    (
        union_coordinates,
        coordinate_to_index,
    ) = build_union_coordinate_space(
        all_coordinates
    )

    n_union = len(union_coordinates)

    print(
        f"Union coordinate count: {n_union}"
    )

    print("Building subject mappings...")

    subject_mappings = (
        build_subject_mappings(
            all_coordinates=all_coordinates,
            coordinate_to_index=coordinate_to_index,
        )
    )

    group_waytotal_sum = None
    group_row_sum = None

    expected_seed_count = None

    for subject_index, (
        connectivity,
        mapping,
        waytotal,
    ) in enumerate(
        zip(
            all_connectivity,
            subject_mappings,
            all_waytotals,
        ),
        start=1,
    ):
        print(
            f"Projecting subject "
            f"{subject_index}/{n_subjects}"
        )

        (
            waytotal_matrix,
            row_matrix,
        ) = build_subject_matrices(
            connectivity=connectivity,
            subject_to_union=mapping,
            n_union=n_union,
            waytotal=waytotal,
        )

        if expected_seed_count is None:
            expected_seed_count = (
                waytotal_matrix.shape[0]
            )

        elif (
            waytotal_matrix.shape[0]
            != expected_seed_count
        ):
            raise ValueError(
                "Subjects have different numbers of cortical "
                "seed rows. "
                f"Expected {expected_seed_count}, but subject "
                f"{subject_index} has "
                f"{waytotal_matrix.shape[0]}"
            )

        if group_waytotal_sum is None:
            group_waytotal_sum = (
                waytotal_matrix.copy()
            )

            group_row_sum = (
                row_matrix.copy()
            )

        else:
            group_waytotal_sum += (
                waytotal_matrix
            )

            group_row_sum += (
                row_matrix
            )

    group_average_waytotal = (
        group_waytotal_sum
        / n_subjects
    )

    group_average_row = (
        group_row_sum
        / n_subjects
    )

    return (
        group_average_waytotal.tocsr(),
        group_average_row.tocsr(),
        union_coordinates,
        subject_mappings,
    )


# Main 

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create group-average structural connectivity "
            "matrices from HCP tractography outputs."
        )
    )

    parser.add_argument(
        "--subjects",
        required=True,
        help=(
            "Text file containing one HCP subject directory "
            "per line."
        ),
    )

    parser.add_argument(
        "--matrix-folder",
        required=True,
        choices=[
            "mat2",
            "mat2_pd",
        ],
        help=(
            "Tractography directory to process. Use 'mat2' for "
            "the standard matrices or 'mat2_pd' for the "
            "path-distance-corrected matrices."
        ),
    )

    parser.add_argument(
        "--output",
        required=True,
        help=(
            "Output filename prefix. For example: "
            "/path/to/group_average_mat2"
        ),
    )

    parser.add_argument(
        "--skip-missing",
        action="store_true",
        help=(
            "Skip subjects with missing or unreadable files "
            "instead of stopping the complete run."
        ),
    )

    args = parser.parse_args()

    subjects_file = Path(
        args.subjects
    )

    if not subjects_file.exists():
        raise FileNotFoundError(
            f"Subject list does not exist: "
            f"{subjects_file}"
        )

    with subjects_file.open("r") as file:
        subject_dirs = [
            Path(line.strip())
            for line in file
            if line.strip()
        ]

    if not subject_dirs:
        raise ValueError(
            "The subject list contains no subject directories"
        )

    print(
        f"Matrix folder: {args.matrix_folder}"
    )

    print(
        f"Found {len(subject_dirs)} subject directories"
    )

    all_coords_L = []
    all_coords_R = []

    all_connectivity_L = []
    all_connectivity_R = []

    all_waytotals_L = []
    all_waytotals_R = []

    included_subjects = []
    skipped_subjects = []

    for subject_index, subject_dir in enumerate(
        subject_dirs,
        start=1,
    ):
        print()
        print(
            f"Loading subject "
            f"{subject_index}/{len(subject_dirs)}: "
            f"{subject_dir.name}"
        )

        try:
            (
                coords_L,
                coords_R,
                connectivity_L,
                connectivity_R,
                waytotal_L,
                waytotal_R,
            ) = load_subject(
                subject_dir=subject_dir,
                matrix_folder=args.matrix_folder,
            )

        except (
            FileNotFoundError,
            OSError,
            ValueError,
        ) as error:
            if args.skip_missing:
                print(
                    f"Skipping {subject_dir.name}: "
                    f"{error}"
                )

                skipped_subjects.append(
                    (
                        subject_dir.name,
                        str(error),
                    )
                )

                continue

            raise

        all_coords_L.append(
            coords_L
        )

        all_coords_R.append(
            coords_R
        )

        all_connectivity_L.append(
            connectivity_L
        )

        all_connectivity_R.append(
            connectivity_R
        )

        all_waytotals_L.append(
            waytotal_L
        )

        all_waytotals_R.append(
            waytotal_R
        )

        included_subjects.append(
            subject_dir.name
        )

        print(
            f"  Left coordinates: "
            f"{coords_L.shape}"
        )

        print(
            f"  Left connectivity entries: "
            f"{connectivity_L.shape}"
        )

        print(
            f"  Left waytotal: "
            f"{waytotal_L}"
        )

        print(
            f"  Right coordinates: "
            f"{coords_R.shape}"
        )

        print(
            f"  Right connectivity entries: "
            f"{connectivity_R.shape}"
        )

        print(
            f"  Right waytotal: "
            f"{waytotal_R}"
        )

    if not included_subjects:
        raise RuntimeError(
            "No valid subjects were loaded"
        )

    print()
    print(
        f"Processing {len(included_subjects)} "
        f"included subjects"
    )

    print()
    print("Processing left hemisphere...")

    (
        average_waytotal_L,
        average_row_L,
        union_coords_L,
        mappings_L,
    ) = process_hemisphere(
        all_coordinates=all_coords_L,
        all_connectivity=all_connectivity_L,
        all_waytotals=all_waytotals_L,
    )

    print()
    print("Processing right hemisphere...")

    (
        average_waytotal_R,
        average_row_R,
        union_coords_R,
        mappings_R,
    ) = process_hemisphere(
        all_coordinates=all_coords_R,
        all_connectivity=all_connectivity_R,
        all_waytotals=all_waytotals_R,
    )

    output_prefix = Path(
        args.output
    )

    output_prefix.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    print()
    print("Saving outputs...")

    save_npz(
        f"{output_prefix}_wt_L.npz",
        average_waytotal_L,
    )

    save_npz(
        f"{output_prefix}_wt_R.npz",
        average_waytotal_R,
    )

    save_npz(
        f"{output_prefix}_row_L.npz",
        average_row_L,
    )

    save_npz(
        f"{output_prefix}_row_R.npz",
        average_row_R,
    )

    np.save(
        f"{output_prefix}_coords_L.npy",
        union_coords_L,
    )

    np.save(
        f"{output_prefix}_coords_R.npy",
        union_coords_R,
    )

    np.save(
        f"{output_prefix}_mapping_L.npy",
        np.asarray(
            mappings_L,
            dtype=object,
        ),
        allow_pickle=True,
    )

    np.save(
        f"{output_prefix}_mapping_R.npy",
        np.asarray(
            mappings_R,
            dtype=object,
        ),
        allow_pickle=True,
    )

    np.savetxt(
        f"{output_prefix}_included_subjects.txt",
        np.asarray(
            included_subjects,
            dtype=str,
        ),
        fmt="%s",
    )

    if skipped_subjects:
        skipped_path = (
            f"{output_prefix}"
            "_skipped_subjects.txt"
        )

        with open(
            skipped_path,
            "w",
        ) as file:
            file.write(
                "Subject\tReason\n"
            )

            for (
                subject_id,
                reason,
            ) in skipped_subjects:
                clean_reason = (
                    reason.replace(
                        "\n",
                        " ",
                    )
                )

                file.write(
                    f"{subject_id}\t"
                    f"{clean_reason}\n"
                )

    print()
    print("Finished.")
    print(
        f"Included subjects: "
        f"{len(included_subjects)}"
    )
    print(
        f"Skipped subjects: "
        f"{len(skipped_subjects)}"
    )
    print(
        f"Waytotal-normalised left shape: "
        f"{average_waytotal_L.shape}"
    )
    print(
        f"Waytotal-normalised right shape: "
        f"{average_waytotal_R.shape}"
    )
    print(
        f"Row-normalised left shape: "
        f"{average_row_L.shape}"
    )
    print(
        f"Row-normalised right shape: "
        f"{average_row_R.shape}"
    )


if __name__ == "__main__":
    main()