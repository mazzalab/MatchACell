def rearrange_sign(raw_data_filtered, annot_file):

    import pandas as pd

    signature_dict = {
        sheet_name: [
            gene.strip()
            for gene in df.iloc[:, 0].dropna().astype(str)
            if gene.strip() in raw_data_filtered.var_names
        ]
        for sheet_name, df in pd.read_excel(
            annot_file, sheet_name=None, header=None
        ).items()
    }

    return signature_dict


def extract_best_res(verdict_file):
    
    import re

    with open(verdict_file, "r") as f:
        text = f.read()

    match = re.search(r"Recommended resolution\s*:\s*(\S+)", text)

    if match:
        return match.group(1)

    raise ValueError(
        f"'Recommended resolution' not found in {verdict_file}"
    )