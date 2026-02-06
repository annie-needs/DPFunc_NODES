# Author: Annie Needs
# Date: January 26th, 2026
# Ojective: a scrip to get the AF2DB .pdb files using the uniprot ids from the ID mapping txt file and map files

import os
import json

with open('../all_id_map.json', 'r') as f:
    data = json.load(f)
i = 1
for key, val in data.items():
    file = f'AF-{val}-F1-model_v6.pdb.gz'
    
    print(i)
    if not os.path.isfile(file):
        source = f'https://alphafold.ebi.ac.uk/files/AF-{val}-F1-model_v6.pdb'
        cmd = f'wget -qO- {source} | gzip > AF-{val}-F1-model_v6.pdb.gz'
        os.system(cmd)
    
    i += 1