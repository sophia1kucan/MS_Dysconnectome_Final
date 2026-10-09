

# Averaging Connectivity and Mapping Dysconnectivity from MS Lesions

**Author:** Sophia Kucan  
**Date:** July 2026

## Overview

This repository contains a computational neuroimaging pipeline for investigating structural dysconnectivity associated with multiple sclerosis (MS) lesions.

The pipeline constructs group-average structural connectivity matrices from diffusion tractography data, maps the lesion-associated dysconnectivity onto the cortical surface, and then quantifies the dysconnectivity using several cortical parcellation atlases.


## Aim

- Construct group-average structural connectivity matrices across subjects.
- Estimate lesion-associated structural dysconnectivity.
- Project dysconnectivity measurements onto the cortical surface.
- Quantify regional dysconnectivity using Yeo-7, Yeo-17 and Glasser parcellations.
- Identify cortical regions exceeding a specified dysconnectivity threshold.
- Construct subject-by-parcel matrices for further statistical analysis.

## Workflow

1. **Group connectivity averaging:** `average_connectivity3.py`
   - Construct group-average connectivity matrices from subject-level tractography data.

2. **Lesion dysconnectivity mapping:** `lesion_dysconnectivity2.py`
   - Combine the lesion masks with group connectivity matrices to generate cortical dysconnectome maps.

3. **Cortical parcellation:** `run_all_parcellations2.sh`
   - Apply Yeo 7, Yeo 17, and Glasser cortical atlases to the dysconnectome maps.

4. **Regional analysis:** `analyse_parcels.py`
   - Identify cortical parcels exceeding a specified dysconnectivity threshold.

5. **Subject-level matrix construction:** `build_dysconnectivity.py`
   - Combine individual parcellated dysconnectomes into a subject-by-parcel matrix for any later analysis.






## Necessary Packages:
To run the scripts, a python virtual environment was used. The necessary packages can be installed from `requirements.txt` using:
```bash
pip install -r requirements.txt
```
### Additional Dependencies:
- Cortical parcellation requires **Connectome Workbench**, which can be accessed via the 'wb_command' command line function.
- The Yeo 7, Yeo 17, and Glasser atals files also need to be available locally.
- File paths in 'run_all_parcellations2.sh' are currently configured for the orginal research environment and will need to be updated accordingly before replicating the pipeline elsewhere.



## Method Summary:
### Group averaging
For each hemisphere:
1. The script creates a union coordinate space of the subjects:
2. The connectivity matrix for each subject is projected onto the union space.
3. The connectivity matrices are averaged:

$$
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


## Regional Dysconnectivity Analysis

Now that individual cortical dysconnectome maps have been generated, the pipeline quantifies dysconnectivity at across different parcellation maps with varying levels of detail.

Three cortical parcellation atlases are used:

- **Yeo 7:** Seven functional networks.
- **Yeo 17:** Seventeen functional networks.
- **Glasser:** 360 cortical parcels.

This means that lesion-associated dysconnectivity can be investigated at different levels of detail.

The following scripts automate parcellation, identify affected regions, and rewrite the results into matrix form for potential statistical analysis.



## Script 3 - Cortical Parcellation
### Aim
'run_all_parcellations2.sh' automates the cortical parcellation and regional dysconnectivity analysis across all subjects.

### The Script:
1. Identifies individual subject dysconnectome files ('.dscalar.nii').
2. Uses Connectome Workbench to parcellate each dysconnecotme using the Yeo 7, Yeo 17, and Glasser atlases.
3. Generates a parcel-level dysconnectivity map ('.pscalar.nii').
4. Executes 'analyse_parcels.py' and uses a dysconnectivity threshold of 0.1.
5. Saves atlas specific affected region reports for each subject.

### Inputs Required
- Subject-level dysconnectome maps (`.dscalar.nii`).
- Yeo 7, Yeo 17 and Glasser atlas files (`.dlabel.nii`).
- Connectome Workbench (`wb_command`).
- `analyse_parcels.py`.

### Outputs
For each subject, the script generates:

- Three parcellated CIFTI files (`.pscalar.nii`).
- Three affected-region reports (`.txt`).




## Script 4 - Identification of Affected Cortical Regions

### Aim
`analyse_parcels.py` identifies cortical regions with dysconnectivity above a defined threshold.

### The script:
1. Loads parcel-level dysconnectivity values from a CIFTI file.
2. Extracts the region labels from the corresponding cortical atlas.
3. Identifies parcels meeting or exceeding the dysconnectivity threshold.
4. Exports the selected parcel IDs, names, and dysconnectivity values.

### Inputs Required
- Parcellated CIFTI file (`.pscalar.nii`).
- Corresponding atlas (`.dlabel.nii`).
- Dysconnectivity threshold.

| Flag | Required | Default | Description |
|---|---|---|---|
| `--map` | Yes | — | Parcellated CIFTI file |
| `--atlas` | Yes | — | Cortical atlas label file |
| `--threshold` | No | 0.0 | Minimum dysconnectivity value |
| `--output` | Yes | — | Output text filename |

### Output
A text file with the parcel ID, region name, and dysconnectivity value for each parcel meeting the threshold.

---

## Script 5 - Subject-by-Parcel Matrix Construction

### Aim
`build_dysconnectivity.py` combines individual subject-level parcellated dysconnectomes into a single matrix for later statistical analysis.

### The script:
1. Loads all matching parcellated CIFTI files from a specified directory.
2. Extracts subject IDs and parcel-level dysconnectivity values.
3. Verifies consistent parcel ordering across subjects.
4. Constructs a matrix with subjects as rows and cortical parcels as columns.
5. Performs quality-control checks.
6. Exports the resulting matrix as a CSV file.

### Inputs Required
- Directory containing parcellated CIFTI files (`.pscalar.nii`).
- Filename prefix identifying the parcellation atlas.

| Flag | Required | Default | Description |
|---|---|---|---|
| `--input-dir` | Yes | — | Directory containing parcellated files |
| `--prefix` | No | `Glasser_` | Filename prefix before subject ID |
| `--output` | Yes | — | Output CSV filename |

### Output
A CSV file containing a subject-by-parcel dysconnectivity matrix.

Each row represents an individual subject, while each column represents a cortical parcel. The entries contain regional dysconnectivity values, providing a structured dataset for any further statistical modelling or investigation of clinical outcomes.






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
- Parcellation is performed separately using the Yeo 7, Yeo 17, and Glasser atlases.
- A dysconnectivity threshold of 0.1 is used to filter affected regions in the processing script.
- The subject by parcel matrix is then used to investigate any association between structural dysconnectivity and outcome.

