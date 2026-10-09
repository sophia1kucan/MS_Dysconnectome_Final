#!/bin/bash

# Input/output folder

INPUT_DIR=/home/pmysk14/ms.lesion.proj/lesion_matrices/final/mat2_row_dysconnectomes


# Parcellation atlases

PARC7=/home/data/MS_lesion_conn/parcellations/Yeo2011_7Networks_N1000.10k_fs_LR.dlabel.nii

PARC17=/home/data/MS_lesion_conn/parcellations/Yeo2011_17Networks_N1000.10k_fs_LR.dlabel.nii

GLASSER=/home/data/MS_lesion_conn/parcellations/Glasser2016.10k_fs_LR.dlabel.nii


# Parcel-analysis Python script

ANALYSE_SCRIPT=/home/pmysk14/ms.lesion.proj/lesion_matrices/analysing_parcels/analyse_parcels.py


# Loop through all dysconnectomes

for FILE in ${INPUT_DIR}/dysconnectome*.dscalar.nii

do

    NAME=$(basename ${FILE})

    ID=${NAME#dysconnectome}
    ID=${ID%.dscalar.nii}


    echo "================================"
    echo "Processing subject ${ID}"
    echo "================================"


    # Yeo 7 networks

    echo "Running Yeo 7 network parcellation..."

    wb_command -cifti-parcellate \
    ${FILE} \
    ${PARC7} \
    COLUMN \
    ${INPUT_DIR}/parc7_${ID}.pscalar.nii


    python ${ANALYSE_SCRIPT} \
    --map ${INPUT_DIR}/parc7_${ID}.pscalar.nii \
    --atlas ${PARC7} \
    --threshold 0.1 \
    --output ${INPUT_DIR}/${ID}_parc7_affected_regions.txt


    # Yeo 17 networks

    echo "Running Yeo 17 network parcellation..."

    wb_command -cifti-parcellate \
    ${FILE} \
    ${PARC17} \
    COLUMN \
    ${INPUT_DIR}/parc17_${ID}.pscalar.nii


    python ${ANALYSE_SCRIPT} \
    --map ${INPUT_DIR}/parc17_${ID}.pscalar.nii \
    --atlas ${PARC17} \
    --threshold 0.1 \
    --output ${INPUT_DIR}/${ID}_parc17_affected_regions.txt


    # Glasser

    echo "Running Glasser parcellation..."

    wb_command -cifti-parcellate \
    ${FILE} \
    ${GLASSER} \
    COLUMN \
    ${INPUT_DIR}/Glasser_${ID}.pscalar.nii


    python ${ANALYSE_SCRIPT} \
    --map ${INPUT_DIR}/Glasser_${ID}.pscalar.nii \
    --atlas ${GLASSER} \
    --threshold 0.1 \
    --output ${INPUT_DIR}/${ID}_Glasser_affected_regions.txt


    echo "Finished ${ID}"

done


echo "================================"
echo "All parcellations finished <3"
echo "================================"