# Process_data.ipynb but in  .py form so I can submit through step 3 as an sbatch job (bc it runs too long and kept getting inturrupted) 

import pandas as pd
import numpy as np
import os
from tqdm.auto import tqdm
#
#
import pickle as pkl

def read_pkl(file_path):
    with open(file_path,'rb') as fr:
        return pkl.load(fr)

def save_pkl(file_path, val):
    fw = open(file_path, 'wb')
    pkl.dump(val, fw)
    fw.close()
#
#
tags = ['mf', 'cc', 'bp']
types = ['train', 'valid', 'test']
#
#
#first, generate pid_list_file
os.makedirs("./processed_file", exist_ok=True)

for tag in tags:
    for tp in types:
        pid_list = set()
        with open(f"./data_dpfunc/{tag}_{tp}_pid_list.txt", 'r') as f:
            lines = f.readlines()
            for line in lines:
                content = line.strip('\n').strip()
                pid_list.add(content)
        pid_list = list(pid_list)
        save_pkl(f"./processed_file/{tag}_{tp}_used_pid_list.pkl", pid_list)
#
#
# Step 2. pid_go_file
for tag in tags:
    for tp in types:
        cmd = f"cp ./data_dpfunc/{tag}_{tp}_go.txt ./processed_file/{tag}_{tp}_go.txt"
        os.system(cmd)
#
#
# pid_pdb_file
from tqdm.auto import tqdm
import gzip
#
#
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

def extract_sequence_and_ca_coords(pdb_file, chain_id=None):
    parser = PDBParser(QUIET=True)
    if pdb_file.endswith('.gz'):
        with gzip.open(pdb_file, 'rt') as gz_file:
            temp_file = pdb_file.replace('.gz', '_temp')
            try:
                with open(temp_file, 'w') as temp:
                    temp.write(gz_file.read())
                
                structure = parser.get_structure('protein', temp_file)
            finally:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
    else:
        structure = parser.get_structure('protein', pdb_file)
    
    results = {}
    
    for model in structure:
        for chain in model:
            if chain_id is None or chain.id == chain_id:
                sequence = ""
                ca_coords = []
                
                for residue in chain:
                    if residue.id[0] == ' ':
                        try:
                            aa = seq1(residue.resname)
                            sequence += aa
                            
                            if 'CA' in residue:
                                ca_atom = residue['CA']
                                coord = ca_atom.get_coord()
                                ca_coords.append((float(coord[0]), float(coord[1]), float(coord[2])))
                            else:
                                print(f"Warning: No CA atom found in residue {residue.resname}{residue.id[1]} of chain {chain.id}")
                                ca_coords.append(None)
                                
                        except KeyError:
                            print(f'Non natural residue in {pdb_file}')
                
                results[chain.id] = {
                    'sequence': sequence,
                    'ca_coords': ca_coords
                }
    
    return results

def extract_single_chain(pdb_file, chain_id='A'):
    results = extract_sequence_and_ca_coords(pdb_file, chain_id)
    
    if chain_id in results:
        return results[chain_id]['sequence'], results[chain_id]['ca_coords']
    else:
        print(f"Chain {chain_id} not found in PDB file")
        return "", []
#
#
pid_list = set()
for tag in tags:
    for tp in types:
        tp_pid_list = read_pkl(f'./processed_file/{tag}_{tp}_used_pid_list.pkl')
        pid_list = pid_list|set(tp_pid_list)
pid_list = list(pid_list)
print(pid_list[0:10])
print(len(pid_list))
#
#
train_id_map = read_pkl('./data_dpfunc/train_id_map.pkl')
valid_id_map = read_pkl('./data_dpfunc/valid_id_map.pkl')
test_id_map = read_pkl('./data_dpfunc/test_id_map.pkl')
#
#
assert len(train_id_map)+len(valid_id_map)+len(test_id_map) == len(set(train_id_map.keys())|set(valid_id_map.keys())|set(test_id_map.keys()))
#
#
all_id_map = {}
for k,v in train_id_map.items():
    all_id_map[k] = v
for k,v in valid_id_map.items():
    all_id_map[k] = v
for k,v in test_id_map.items():
    all_id_map[k] = v
#
#
for k,v in list(all_id_map.items())[:10]:
    print(k, v)
print(len(all_id_map))
#
#
tp_map = pd.read_table('./idmapping_2025_06_13.tsv')
tp_map.head()
#
#
for idx, row in tp_map.iterrows():
    all_id_map[row['Entry Name']] = row['Entry']
#
#
import json

with open('all_id_map.json', 'w') as json_file:
    json.dump(all_id_map, json_file, indent=4)
print("dictionary has been written to all_id_map.json")
#
#
pdb_points_info = {}
pdb_seq_info = {}
unseen_proteins = set()

for protein in tqdm(pid_list):
    uni_id = all_id_map[protein]
    pdb_file = f"./AF2DB/AF-{uni_id}-F1-model_v6.pdb.gz"
    if not os.path.exists(pdb_file):
        unseen_proteins.add(protein)
        continue
    
    sequence, coords_list = extract_single_chain(pdb_file, 'A')
    
    if sequence and coords_list:
        valid_coords = [coord for coord in coords_list if coord is not None]
        valid_sequence = ''.join([sequence[i] for i, coord in enumerate(coords_list) if coord is not None])
        
        if len(valid_coords)==0 or len(valid_sequence)==0:
            print(f"Empty {protein}, {uni_id}")
        
        pdb_points_info[protein] = valid_coords
        pdb_seq_info[protein] = valid_sequence
#
#
save_pkl('./processed_file/pdb_points.pkl', pdb_points_info)
save_pkl('./processed_file/pdb_seqs.pkl', pdb_seq_info)
save_pkl('./processed_file/unseen_proteins.pkl', unseen_proteins)
#
#
# DO step 4 later bc process_esm.py needs other files ran. 