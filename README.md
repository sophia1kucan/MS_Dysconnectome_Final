# Averaging Connectivity and Mapping Dysconnectivity from Lesions
## Author: Sophia Kucan
## July 2026

This README describes the purpose and how to use the scripts `average_connectivity.py` and `lesion_dysconnectivity.py` from the folder `final`.

## Aim:

- Compute average connectivity across subjects or regions.
- Map dysconnectivity patterns arising from lesions.
- Generate readable outputs for downstream analysis.

## Workflow:

1. Subject connectivity data.
2. Use `average_connectivity3.py`.
3. Group average connectivity matrices.
4. Use `lesion_dysconnectivity2.py`.
5. CIFTI lesion dysconnectome.

## Necessary Packages:
To run the scripts, a python virtual environment was used. The necessary packages can be installed from `requirements.txt` using:
```bash
pip install -r requirements.txt
```


## Method Summary:
### Group averaging
For each hemisphere:
1. The script creates a union coordinate space of the subjects:
2. The connectivity matrix for each subject is projected onto the union space.
3. The connectivity matrices are averaged:

\text{Group average connectivity} =
\frac{\sum_{i=1}^{N}\text{Connectivity matrix}_{i}}
{N}
$$


### Lesion dysconnectivity
The lesion dysconnectivity stage takes the average connectivity matrix from the previous section and combines this with a lesion mask to identify the loss of connectivity caused by the lesion. This is then projected back onto the cortical surface.
### General pipeline:
1. Lesion mask.
2. Resample into tract space.
3. Tract space lookup.
4. Finding the connectivity associated with the lesion.
5. Calculate the average dysconnectivity caused by the lesion.
6. Project this back onto the cortical surface.
7. CIFTI map.

## Script 1 - Group Connectivity Averaging
### Aim:
`average_connectivity3.py` calculates the average connectivity matrix of a group of subjects.
The script:
1. Loads the connectivity data.
2. Builds a union coordinate space across all subjects.
3. Projects each subject's connectivity matrix onto the union space.
4. Takes the average of the connectivity matrix across subjects.
5. Saves the resulting group matrices and mapping information.

### Input Data Requirements:
Each subject directory should contain:
`subject_directory`
- `omatrix2`
- `omatrix2_L.white.resampled`

The script expects:
- Left hemisphere coordinate file.
- Right hemisphere coordinate file.
- Sparse tractography connectivity matrices.

| Flag | Required | Default | Description |
|---|---|---|---|
| `--subjects` | Yes | — | Path to a text file listing subject directories, one per line |
| `--output` | Yes | — | Prefix used for all output files (e.g., `group_average`) |


### How to Run the Averaging Script:
You will need a text file that contains subject directories, where each line represents one subject. 

#### Example: 

```bash
subjects.txt
``` 

#### Contents:

```bash
/home/data/MS_lesion_conn/HCP_example_subjects/******
/home/data/MS_lesion_conn/HCP_example_subjects/******
/home/data/MS_lesion_conn/HCP_example_subjects/******
/home/data/MS_lesion_conn/HCP_example_subjects/******
/home/data/MS_lesion_conn/HCP_example_subjects/******
/home/data/MS_lesion_conn/HCP_example_subjects/******
/home/data/MS_lesion_conn/HCP_example_subjects/******
/home/data/MS_lesion_conn/HCP_example_subjects/******
/home/data/MS_lesion_conn/HCP_example_subjects/******
/home/data/MS_lesion_conn/HCP_example_subjects/******
```

Each line represents one subject.

#### Run:
```bash
python average_connectivity.py
--subjects subjects.txt
--output group_average
```

### Outputs
The script generates the files:

```bash
group_average_L.npz
group_average_R.npz
```
These are the average connectivity matrices of the sparse group matrices.

```bash
group_average_coords_L.npy
group_average_coords_R.npy
```
These are the coordinates that index the shared tract space.

```bash
group_average_mapping_L.npy
group_average_mapping_R.npy
```
These are the mappings which are needed in the later lesion dysconnectivity calculations.

## Script 2 - Lesion Dysconnectivity Mapping

### Aim:
`lesion_dysconnectivity2.py` calculates the lack of connectivity caused by a lesion in a given area of the brain.

### The script:
1. Loads the group average connectivity matrices computed in `average_connectivity.py`.
2. Loads the lesion mask.
3. Converts the voxels affected by the lesion into tract space coordinates.
4. Removes the connectivity values affected by the lesion.
5. Projects the resulting matrix back onto the cortical surface.
6. Saves this as a CIFTI dysconnectome map.

### Inputs Required
The script requires:
- Group average connectivity matrices (generated in the first script):
`group_average_L.npz` \
`group_average_R.npz` 
- Subjects mappings (generated in the first script):
`group_average_mapping_L.npy` \
`group_average_mapping_R.npy` 
- Lesion mask:
A lesion mask `lesion_mask.nii.gz`. 
- Reference subject:
A reference subject that contains:
  - Cortical surfaces
  - Atlas masks
  - Tract space lookup images

| Flag | Required | Default | Description |
|---|---|---|---|
| `--matrix-L` | Yes | — | Path to the left-hemisphere group average matrix (`group_average_L.npz`) |
| `--matrix-R` | Yes | — | Path to the right-hemisphere group average matrix (`group_average_R.npz`) |
| `--mapping-L` | Yes | — | Path to the left-hemisphere group average mapping (`group_average_mapping_L.npy`) |
| `--mapping-R` | Yes | — | Path to the right-hemisphere group average mapping (`group_average_mapping_R.npy`) |
| `--reference-subject` | Yes | — | Path to a reference subject directory containing cortical surfaces, atlas masks, and tract space lookup images |
| `--lesion` | Yes | — | Path to the lesion mask (`lesion_mask.nii.gz`) |
| `--output` | Yes | — | Output filename for the resulting CIFTI dysconnectome map (e.g., `lesion_disconnectome.dscalar.nii`) |

### How to Run the Lesion Mapping Script
#### Run:

```bash
python lesion_dysconnectivity2.py \
--matrix-L group_average_L.npz \
--matrix-R group_average_R.npz \
--mapping-L group_average_mapping_L.npy \
--mapping-R group_average_mapping_R.npy \
--reference-subject /path/to/reference_subject \
--lesion /path/to/lesion_mask.nii.gz \
--output lesion_disconnectome.dscalar.nii
```

### Output:
The final output is: \
`lesion_disconnectome.dscalar.nii`\
Which can then be viewed using Connectome Workbench.
Open the CIFTI file:
`wb_view lesion_disconnectome.dscalar.nii`.
I would recommend viewing this on top of inflated brain maps.

## Troubleshooting
### File not found errors
Check that:
- The required files exist.
- The required subject directories are correct.
- The reference subject contains the required lookup files.

### Lesion does not overlap tract space
If the script shows:\
`Lesion does not overlap any tract space coordinates`.\
Check:
- Orientation of the lesion image.
- The reference lookup image.
- The lesion registration.

## Notes
Connectivity averaging is done separately from lesion maps to allow:
- Average group connectivity matrices to be reused.
- Multiple lesions to be processed without needing to recompute averages.
- Broader application amongst different subject groups.

