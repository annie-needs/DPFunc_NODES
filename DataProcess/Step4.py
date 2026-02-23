# Python script to run step 4 of data processing


import pandas as pd
import numpy as np
import os
from tqdm.auto import tqdm
import pickle as pkl
import dgl
import math
import torch

def read_pkl(file_path):
    with open(file_path,'rb') as fr:
        return pkl.load(fr)

def save_pkl(file_path, val):
    fw = open(file_path, 'wb')
    pkl.dump(val, fw)
    fw.close()

def get_dis(point1, point2):
    dis_x = point1[0] - point2[0]
    dis_y = point1[1] - point2[1]
    dis_z = point1[2] - point2[2]
    return math.sqrt(dis_x*dis_x + dis_y*dis_y + dis_z*dis_z)

def process_input_pdb_file(tag, part, pid_list, pdb_points_info, pdb_seqs, thresholds=12):
    protein_map = read_pkl('./processed_file/protein_map.pkl')
    pdb_graphs = []
    p_cnt = 0
    file_idx = 0
    for pid in tqdm(pid_list):
        p_cnt += 1
        points = pdb_points_info[pid]
        
        u_list = []
        v_list = []
        dis_list = []
        for uid, amino_1 in enumerate(points):
            for vid, amino_2 in enumerate(points):
                if uid==vid:
                    continue
                dist = get_dis(amino_1, amino_2)
                if dist<=thresholds:
                    u_list.append(uid)
                    v_list.append(vid)
                    dis_list.append(dist)
        u_list, v_list = torch.tensor(u_list), torch.tensor(v_list)
        dis_list = torch.tensor(dis_list)

        graph = dgl.graph((u_list, v_list), num_nodes=len(points))
        graph.edata['dis'] = dis_list

        # graph node feature - esm
        esm_file_idx = protein_map[pid]
        esm_features = read_pkl(f"./processed_file/esm_emds/esm_part_{esm_file_idx}.pkl")
        node_features = esm_features[pid]
        assert node_features.shape[0]==graph.num_nodes()
        graph.ndata['x'] = torch.from_numpy(node_features)
        pdb_graphs.append(graph)

        if p_cnt%5000==0:
            save_pkl('./processed_file/graph_features/{}_{}_whole_pdb_part{}.pkl'.format(tag, part, file_idx), pdb_graphs)
            p_cnt = 0
            file_idx += 1
            pdb_graphs = []
    if len(pdb_graphs)>0:
        save_pkl('./processed_file/graph_features/{}_{}_whole_pdb_part{}.pkl'.format(tag, part, file_idx), pdb_graphs)
    return file_idx

tags = ['mf', 'cc', 'bp']
types = ['train', 'valid', 'test']

os.makedirs('./processed_file/graph_features/', exist_ok=True)

all_protein_list = []
for tag in tags:
    for tp in types:
        pid_list = read_pkl(f"./processed_file/{tag}_{tp}_used_pid_list.pkl")
        all_protein_list+=pid_list
all_protein_list = list(set(all_protein_list))

pdb_points_info = read_pkl('./processed_file/pdb_points.pkl')
pdb_seqs = read_pkl('./processed_file/pdb_seqs.pkl')


to_remove = [pid for pid in all_protein_list if pid not in pdb_points_info or pid not in pdb_seqs]

print(len(to_remove))
print(len(all_protein_list))

all_protein_list = [pid for pid in all_protein_list if pid not in to_remove]
        
for tag in tags:
    for tp in types:
        residual_protein_list = []
        pid_list = read_pkl(f"./processed_file/{tag}_{tp}_used_pid_list.pkl")
        print(f'inital count = {len(pid_list)}')

        residual_protein_list = [pid for pid in pid_list if pid not in to_remove]
        print(f'final count = {len(residual_protein_list)}')

        save_pkl(f'./processed_file/{tag}_{tp}_real_used_pid_list.pkl', residual_protein_list)

print(len(all_protein_list))


assert len(set(pdb_points_info.keys()))==len(set(pdb_seqs.keys()))
assert len(set(pdb_points_info.keys())&set(all_protein_list))==len(all_protein_list)
assert len(set(pdb_seqs.keys())&set(all_protein_list))==len(all_protein_list)

for tag in tags:
    if tag=='mf':
        continue
    for tp in types:
        pid_list = read_pkl(f"./processed_file/{tag}_{tp}_real_used_pid_list.pkl")
        max_cnt = process_input_pdb_file(tag, tp, pid_list, pdb_points_info, pdb_seqs)
        if tp=='train':
            print(f"{tag}-{tp}-train_file_count-{max_cnt}")