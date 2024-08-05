import re
import pandas as pd


def extract_pca_contributions(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Extract properties
    properties_line = lines[12]
    properties = eval(properties_line.strip())

    # Extract PCA contributions starting lines
    pca1_start = 62
    pca2_start = 105
    pca2_end = 146
    # Extract PCA1 contributions
    pca1_contributions = {}
    for line in lines[pca1_start + 1 : pca2_start]:
        match = re.search(
            r"([\+\-]?\d+\.\d+e[\+\-]\d+).*?/\s+([\+\-]?\d+\.\d+e[\+\-]\d+)\s+\*\s+\((.*?)\s+-",
            line,
        )
        if match:
            contribution = float(match.group(1))
            z_score = float(match.group(2))
            property_name = match.group(3)
            pca1_contributions[property_name] = contribution

    # Extract PCA2 contributions
    pca2_contributions = {}
    for line in lines[pca2_start + 1 : pca2_end]:
        match = re.search(
            r"([\+\-]?\d+\.\d+e[\+\-]\d+).*?/\s+([\+\-]?\d+\.\d+e[\+\-]\d+)\s+\*\s+\((.*?)\s+-",
            line,
        )
        if match:
            contribution = float(match.group(1))
            z_score = float(match.group(2))
            property_name = match.group(3)
            pca2_contributions[property_name] = contribution

    # Create DataFrame
    data = {
        "Property": properties,
        "PCA1contribution": [pca1_contributions.get(prop, 0) for prop in properties],
        "PCA2contribution": [pca2_contributions.get(prop, 0) for prop in properties],
    }
    df = pd.DataFrame(data)
    return df


df_contributions = extract_pca_contributions("../data/pca_data.txt")
df_contributions.to_csv(
    "../data/pca_contributions.csv", index=False, float_format="%.10g"
)
