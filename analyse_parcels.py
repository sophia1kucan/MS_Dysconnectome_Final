import numpy as np
import nibabel as nib
import argparse

def load_parcel_values(pscalar_file):
    img = nib.load(pscalar_file)
    data = img.get_fdata()
    values = data.squeeze()
    return values


# def load_labels(label_file):
#     labels = {}

#     with open(label_file) as f:
#         lines = [line.strip() for line in f if line.strip()]
#     for i in range(0, len(lines), 2):
#         region_name = lines[i]
#         id_line = lines[i + 1]
#         parcel_id = int(id_line.split()[0])
#         labels[parcel_id] = region_name
#     return labels

def load_labels_from_dlabel(dlabel_file):
    img = nib.load(dlabel_file)

    label_axis = img.header.get_axis(0)

    labels = {}

    for parcel_id, label_info in label_axis.label[0].items():
        if parcel_id == 0:
            continue

        labels[parcel_id] = label_info[0]

    return labels

def find_affected_regions(values, labels, threshold):
    affected = []
    for i, value in enumerate(values):
        parcel_id = i + 1
        if value >= threshold:
            region = labels.get(
                parcel_id,
                "Unknown"
            )

            affected.append({
                "parcel": parcel_id,
                "region":region,
                "disconnectivity": value
            })
    return affected

def main():
    parser = argparse.ArgumentParser(
        description="Identify affected cortical parcels from disconnectome"
    )

    # parser.add_argument(
    #     "--map",
    #     required=True,
    #     help="Glasser labels txt file"
    # )

    # parser.add_argument(
    #     "--labels",
    #     required=True,
    #     help="Glasser labesl txt file"

    # )

    parser.add_argument(
        "--map",
        required=True,
        help="Parcellated CIFTI pscalar file"
    )

    parser.add_argument(
        "--atlas",
        required=True,
        help="CIFTI dlabel atals file"
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=0.0,
        help="Dysconncetivity threshold"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output text file"
    )

    args = parser.parse_args()

    values = load_parcel_values(args.map)
    labels = load_labels_from_dlabel(args.atlas)

    affected = find_affected_regions(
        values,
        labels,
        args.threshold
    )

    with open(args.output, "w") as f:
        f.write(
            "Parcel\tRegion\tDisconnectivity\n"
        )

        for region in affected:
            f.write(
                f"{region['parcel']}\t"
                f"{region['region']}\t"
                f"{region['disconnectivity']:.5f}\n"
            )
    
    print(f"Found {len(affected)} affected parcels")

if __name__ == "__main__":
    main()