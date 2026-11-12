# Calculate actifpTM for AlphaFold3 predictions
# Runs reactifptm on every AlphaFold3 output folder to calculate actifpTM, an interface score
# restricted to confidently predicted residues (Varga, Ovchinnikov & Schueler-Furman, 2025).
#
# Input:
#   - AF3_OUTPUT: folder of AlphaFold3 outputs (from run_alphafold3.sh), one subfolder per protein pair,
#                 each containing <name>_confidences.json and <name>_model.cif
#
# Output:
#   - results.json in each subfolder, with the actifpTM scores for that pair


import subprocess
from pathlib import Path

# path to reactifptm; change to the full path if it is not on the PATH
REACTIFPTM = "reactifptm"

# AlphaFold3 output folder
AF3_OUTPUT = Path("af3_output")

for folder in AF3_OUTPUT.iterdir():
    if not folder.is_dir():
        continue
    name = folder.name
    subprocess.run([
        REACTIFPTM,
        folder / f"{name}_confidences.json",
        folder / f"{name}_model.cif",
        "-o", folder / "results.json",
    ])
