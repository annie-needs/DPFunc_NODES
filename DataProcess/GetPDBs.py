# Author: Annie Needs
# Date: January 26th, 2026
# Ojective: a scrip to get the AF2DB .pdb files using the uniprot ids from the ID mapping txt file and map files

import os
import json

with open('../all_id_map.json', 'r') as f:
    data = json.load(f)

i=1

for key, val in data.items():
    source = f'https://alphafold.ebi.ac.uk/files/AF-{val}-F1-model_v6.pdb'
    if i > 10:
        break
    i += 1

    cmd = f'wget {source}'
    os.system(cmd)

