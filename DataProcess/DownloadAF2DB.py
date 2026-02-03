# dowloading AF2DB pdb files
# used ChatGPT to help write script

import pickle
from pathlib import Path

pkl_files = [
    './data_dpfunc/train_id_map.pkl',
    './data_dpfunc/test_id_map.pkl',
    './data_dpfunc/valid_id_map.pkl'
]

uniprot_ids = set()

for pkl in pkl_files:
    with open (pkl,'rb') as f:
        data = pickle.load(f) 

    if isinstance(data,dict):
        for v in data.values():
            if isinstance(v, str):
                uniprot_ids.add(v.strip())

    elif isinstance(data, list):
        for item in data:
            if isinstance(item, dict):
                for k in ("Entry","UniProt","uniprot_id"):
                    if k in item:
                        uniprot_ids.add(str(item[k]).strip())

print(f'collected {len(uniprot_ids)} unique uniprot ids')

def normalize_ids(uid):
    uid = uid.strip()
    return uid.split('-')[0]

uniprot_ids = {normalize_ids(uid) for uid in uniprot_ids}

import requests
BASE = 'https://alphafold.ebi.ac.uk/files'

def afdb_exists(uniprot):
    af_id = f'AF-{uniprot}-F1'
    url = f'{BASE}/{af_id}-model_v4.pdb.gz'
    
    r = requests.head(url,timeout=10)
    return r.status_code ==200

valid =set()
missing = set()

for u in uniprot_ids:
    if afdb_exists(u):
        valid.add(u)
    else:
        missing.add(u)
    
print(f'afdb available: {len(valid)}')
print(f'afdb missing: {len(missing)}')

