import argparse
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.cifti2 import (
    BrainModelAxis,
    Cifti2Image,
    ScalarAxis,
)
from nibabel.processing import resample_from_to
from scipy.sparse import load_npz


#Load the reference data

def load_reference_data(
    surface_reference,
    tractography_reference,
    matrix_folder,
):


    surface_reference = Path(surface_reference)
    tractography_reference = Path(tractography_reference)

  

    atlas_L_path = (
        surface_reference
        / "atlasroi.resampled.L.shape.gii"
    )

    atlas_R_path = (
        surface_reference
        / "atlasroi.resampled.R.shape.gii"
    )

  

    surface_L_path = (
        surface_reference
        / "L.white.resampled.surf.gii"
    )

    surface_R_path = (
        surface_reference
        / "R.white.resampled.surf.gii"
    )


    lookup_L_path = (
        tractography_reference
        / matrix_folder
        / "omatrix2_L.white.resampled"
        / "lookup_tractspace_fdt_matrix2.nii.gz"
    )

    lookup_R_path = (
        tractography_reference
        / matrix_folder
        / "omatrix2_R.white.resampled"
        / "lookup_tractspace_fdt_matrix2.nii.gz"
    )

    required_files = [
        atlas_L_path,
        atlas_R_path,
        surface_L_path,
        surface_R_path,
        lookup_L_path,
        lookup_R_path,
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
            "Missing required reference files:\n"
            f"{missing_text}"
        )

    #Load cortical masks
    atlas_L = nib.load(atlas_L_path)
    atlas_R = nib.load(atlas_R_path)

    mask_L = atlas_L.darrays[0].data.astype(bool)
    mask_R = atlas_R.darrays[0].data.astype(bool)

    #Load surfaces for vertex counts
    surface_L = nib.load(surface_L_path)
    surface_R = nib.load(surface_R_path)

    n_vertices_L = surface_L.darrays[0].data.shape[0]
    n_vertices_R = surface_R.darrays[0].data.shape[0]

    #Load tract lookup images
    lookup_L = nib.load(lookup_L_path)
    lookup_R = nib.load(lookup_R_path)

    return (
        mask_L,
        mask_R,
        n_vertices_L,
        n_vertices_R,
        lookup_L,
        lookup_R,
    )


#Map to reference subject

def load_reference_mapping(
    mapping_file,
    included_subjects_file,
    reference_subject_id,
):


    mapping_file = Path(mapping_file)
    included_subjects_file = Path(included_subjects_file)

    if not mapping_file.exists():
        raise FileNotFoundError(
            f"Mapping file does not exist: {mapping_file}"
        )

    if not included_subjects_file.exists():
        raise FileNotFoundError(
            "Included-subject list does not exist: "
            f"{included_subjects_file}"
        )

    mappings = np.load(
        mapping_file,
        allow_pickle=True,
    )

    with included_subjects_file.open("r") as file:
        subject_entries = [
            line.strip()
            for line in file
            if line.strip()
        ]

    subject_ids = [
        Path(entry).name
        for entry in subject_entries
    ]

    if len(mappings) != len(subject_ids):
        raise ValueError(
            "The number of mappings does not match the number "
            "of subjects in the included-subject file. "
            f"Mappings: {len(mappings)}, "
            f"subjects: {len(subject_ids)}"
        )

    reference_subject_id = Path(
        str(reference_subject_id)
    ).name

    if reference_subject_id not in subject_ids:
        raise ValueError(
            f"Reference subject {reference_subject_id} was not "
            f"found in {included_subjects_file}"
        )

    mapping_index = subject_ids.index(
        reference_subject_id
    )

    mapping = np.asarray(
        mappings[mapping_index],
        dtype=np.int64,
    )

    print(
        f"Using mapping for reference subject "
        f"{reference_subject_id}, index {mapping_index}"
    )

    return mapping

  


#Mapping the lesions to tract space 

def lesion_to_union(
    lesion_img,
    lookup_img,
    subject_mapping,
):


    lesion_data = lesion_img.get_fdata()

    # Treat every non-zero voxel as part of the lesion.
    lesion_mask = lesion_data != 0

    lesion_voxels = np.argwhere(
        lesion_mask
    )

    if lesion_voxels.shape[0] == 0:
        raise ValueError(
            "The lesion mask contains no non-zero voxels"
        )

    lookup = lookup_img.get_fdata().astype(
        np.int64
    )

    #Lesion and lookup images are on the same grid because the
    # lesion has already been resampled to the lookup image.
    subject_indices = lookup[
        lesion_voxels[:, 0],
        lesion_voxels[:, 1],
        lesion_voxels[:, 2],
    ]

    # Lookup value zero is background.
    subject_indices = subject_indices[
        subject_indices > 0
    ]

    if subject_indices.size == 0:
        return np.array(
            [],
            dtype=np.int64,
        )

    #FSL indices are one-based
    subject_indices = (
        subject_indices - 1
    )

    #Remove indices outside the reference-subject mapping
    valid_indices = (
        (subject_indices >= 0)
        & (subject_indices < len(subject_mapping))
    )

    subject_indices = subject_indices[
        valid_indices
    ]

    if subject_indices.size == 0:
        return np.array(
            [],
            dtype=np.int64,
        )

    union_indices = subject_mapping[
        subject_indices
    ]

    return np.unique(
        union_indices
    )


def calculate_lesion_map(
    group_average,
    lesion_img,
    lookup_img,
    subject_mapping,
    hemisphere_name,
):


    # Resample the lesion onto the tract space reference image

    lesion_resampled = resample_from_to(
        lesion_img,
        lookup_img,
        order=0,
    )

    lesion_indices = lesion_to_union(
        lesion_img=lesion_resampled,
        lookup_img=lookup_img,
        subject_mapping=subject_mapping,
    )

    if lesion_indices.size == 0:
        print(
            f"  Warning: lesion does not overlap valid "
            f"{hemisphere_name} tract-space locations."
        )

        return np.zeros(
            group_average.shape[0],
            dtype=np.float32,
        )

    if lesion_indices.max() >= group_average.shape[1]:
        raise ValueError(
            f"A {hemisphere_name} lesion union index exceeds "
            "the number of columns in the group-average matrix"
        )

    print(
        f"  {hemisphere_name} overlapping tract-space "
        f"locations: {len(lesion_indices)}"
    )

    lesion_matrix = group_average[
        :,
        lesion_indices,
    ]

    # Average connectivity across all tract-space locations
    # intersected by the lesion.
    lesion_map = np.asarray(
        lesion_matrix.mean(axis=1)
    ).ravel()

    return lesion_map.astype(
        np.float32
    )


#Surface expansion

def expand_to_surface(
    data,
    mask,
    hemisphere_name,
):
    """
    Expand cortical seed values to the full surface, putting
    zeros in medial-wall vertices.
    """

    data = np.asarray(
        data,
        dtype=np.float32,
    )

    expected_values = int(
        mask.sum()
    )

    if len(data) != expected_values:
        raise ValueError(
            f"{hemisphere_name} map length does not match the "
            f"number of non-medial-wall vertices. "
            f"Map values: {len(data)}; "
            f"mask vertices: {expected_values}"
        )

    full_surface = np.zeros(
        len(mask),
        dtype=np.float32,
    )

    full_surface[mask] = data

    return full_surface


#Save as CIFTI

def save_cifti(
    left_map,
    right_map,
    n_vertices_L,
    n_vertices_R,
    output,
    map_name,
):


    output = Path(output)

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    left_map = np.asarray(
        left_map,
        dtype=np.float32,
    )

    right_map = np.asarray(
        right_map,
        dtype=np.float32,
    )

    if len(left_map) != n_vertices_L:
        raise ValueError(
            "Left surface-map length does not match the left "
            "surface vertex count"
        )

    if len(right_map) != n_vertices_R:
        raise ValueError(
            "Right surface-map length does not match the right "
            "surface vertex count"
        )

    data = np.concatenate(
        [
            left_map,
            right_map,
        ]
    )[np.newaxis, :]

    left_axis = BrainModelAxis.from_surface(
        np.arange(n_vertices_L),
        nvertex=n_vertices_L,
        name="CIFTI_STRUCTURE_CORTEX_LEFT",
    )

    right_axis = BrainModelAxis.from_surface(
        np.arange(n_vertices_R),
        nvertex=n_vertices_R,
        name="CIFTI_STRUCTURE_CORTEX_RIGHT",
    )

    brain_axis = (
        left_axis
        + right_axis
    )

    scalar_axis = ScalarAxis(
        [map_name]
    )

    image = Cifti2Image(
        dataobj=data,
        header=(
            scalar_axis,
            brain_axis,
        ),
    )

    image.to_filename(
        str(output)
    )


# Main 

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate a bilateral lesion dysconnectivity CIFTI "
            "map using group-average HCP tractography matrices."
        )
    )

    parser.add_argument(
        "--matrix-L",
        required=True,
        help=(
            "Left hemisphere group-average connectivity matrix "
            "in .npz format."
        ),
    )

    parser.add_argument(
        "--matrix-R",
        required=True,
        help=(
            "Right hemisphere group-average connectivity matrix "
            "in .npz format."
        ),
    )

    parser.add_argument(
        "--mapping-L",
        required=True,
        help=(
            "Left hemisphere subject-to-union mapping file "
            "in .npy format."
        ),
    )

    parser.add_argument(
        "--mapping-R",
        required=True,
        help=(
            "Right hemisphere subject-to-union mapping file "
            "in .npy format."
        ),
    )

    parser.add_argument(
        "--included-subjects",
        required=True,
        help=(
            "Text file listing the HCP subjects included when "
            "the average matrices were generated."
        ),
    )

    parser.add_argument(
        "--surface-reference",
        required=True,
        help=(
            "HCP reference folder containing atlasroi masks and "
            "surface files. For example: "
            "/home/data/MS_lesion_conn/"
            "HCP_example_subjects/100206"
        ),
    )

    parser.add_argument(
        "--tractography-reference",
        required=True,
        help=(
            "HCP reference folder containing mat2 and mat2_pd. "
            "For example: "
            "/home/data/MS_lesion_conn/"
            "HCP_fdt_matrices/100206"
        ),
    )

    parser.add_argument(
        "--reference-subject-id",
        default="100206",
        help=(
            "Subject ID whose lookup images and mapping should "
            "be used. Default: 100206."
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
            "Use 'mat2' for standard tractography or 'mat2_pd' "
            "for path-distance-corrected tractography."
        ),
    )

    parser.add_argument(
        "--lesion",
        required=True,
        help=(
            "Input real or fake lesion mask in .nii or .nii.gz "
            "format."
        ),
    )

    parser.add_argument(
        "--output",
        required=True,
        help=(
            "Output CIFTI filename ending in .dscalar.nii."
        ),
    )

    parser.add_argument(
        "--map-name",
        default="lesion_dysconnectome",
        help=(
            "Name stored inside the output CIFTI map."
        ),
    )

    args = parser.parse_args()

  

    print("Loading connectivity matrices...")

    group_average_L = load_npz(
        args.matrix_L
    )

    group_average_R = load_npz(
        args.matrix_R
    )

    print(
        f"Left matrix shape: "
        f"{group_average_L.shape}"
    )

    print(
        f"Right matrix shape: "
        f"{group_average_R.shape}"
    )



    print("Loading reference-subject mappings...")

    mapping_L = load_reference_mapping(
        mapping_file=args.mapping_L,
        included_subjects_file=args.included_subjects,
        reference_subject_id=args.reference_subject_id,
    )

    mapping_R = load_reference_mapping(
        mapping_file=args.mapping_R,
        included_subjects_file=args.included_subjects,
        reference_subject_id=args.reference_subject_id,
    )

 

    print("Loading reference data...")

    (
        mask_L,
        mask_R,
        n_vertices_L,
        n_vertices_R,
        lookup_L,
        lookup_R,
    ) = load_reference_data(
        surface_reference=args.surface_reference,
        tractography_reference=args.tractography_reference,
        matrix_folder=args.matrix_folder,
    )

    print(
        f"Left surface vertices: {n_vertices_L}"
    )

    print(
        f"Right surface vertices: {n_vertices_R}"
    )

    print(
        f"Left non-medial-wall vertices: {mask_L.sum()}"
    )

    print(
        f"Right non-medial-wall vertices: {mask_R.sum()}"
    )


    print("Loading lesion mask...")

    lesion_path = Path(
        args.lesion
    )

    if not lesion_path.exists():
        raise FileNotFoundError(
            f"Lesion mask does not exist: {lesion_path}"
        )

    lesion_img = nib.load(
        lesion_path
    )

    lesion_data = lesion_img.get_fdata()

    print(
        f"Lesion image shape: {lesion_img.shape}"
    )

    print(
        f"Non-zero lesion voxels before resampling: "
        f"{np.count_nonzero(lesion_data)}"
    )


    print("Calculating left hemisphere dysconnectivity...")

    left_map = calculate_lesion_map(
        group_average=group_average_L,
        lesion_img=lesion_img,
        lookup_img=lookup_L,
        subject_mapping=mapping_L,
        hemisphere_name="left",
    )

    print("Calculating right hemisphere dysconnectivity...")

    right_map = calculate_lesion_map(
        group_average=group_average_R,
        lesion_img=lesion_img,
        lookup_img=lookup_R,
        subject_mapping=mapping_R,
        hemisphere_name="right",
    )

 

    print("Expanding maps to full cortical surfaces...")

    left_full = expand_to_surface(
        data=left_map,
        mask=mask_L,
        hemisphere_name="Left",
    )

    right_full = expand_to_surface(
        data=right_map,
        mask=mask_R,
        hemisphere_name="Right",
    )

    # Save as CIFTI

    print("Saving bilateral CIFTI map...")

    save_cifti(
        left_map=left_full,
        right_map=right_full,
        n_vertices_L=n_vertices_L,
        n_vertices_R=n_vertices_R,
        output=args.output,
        map_name=args.map_name,
    )

    print("Finished.")

    print(
        f"Left range: "
        f"{left_map.min():.8g} to "
        f"{left_map.max():.8g}"
    )

    print(
        f"Right range: "
        f"{right_map.min():.8g} to "
        f"{right_map.max():.8g}"
    )

    print(
        f"Output saved to: "
        f"{Path(args.output).resolve()}"
    )


if __name__ == "__main__":
    main()