import uuid
import glob
import os

import numpy as np


def merge_files(data_folder):

    for subfolder in glob.glob(os.path.join(data_folder, '*', '*')):

        normal_files = glob.glob(os.path.join(subfolder, 'data-subset-0-*.npz'))
        abnormal_files = glob.glob(os.path.join(subfolder, 'data-subset-abnormal-*.npz'))

        merge_subfolder(list_of_files=normal_files, save_folder=subfolder)
        merge_subfolder(list_of_files=abnormal_files, save_folder=subfolder, suffix='abnormal')

def merge_subfolder(list_of_files, save_folder, suffix='0'): 
    
    numberings, idxs, seq_counts = [], [], 0
    
    for data in [np.load(fname, allow_pickle=True) for fname in list_of_files]:
        
        numberings.append(data['numberings'])
        idxs.append(data['idxs'])
        seq_counts += data['idxs'].shape[0]
        
        if seq_counts > 50_000_000:
            
            unique_id = str(uuid.uuid4())
            save_file = os.path.join(save_folder, "data-subset-{}-{}.npz".format(suffix, unique_id))
            np.savez_compressed(save_file, 
                                numberings=np.concatenate(numberings), 
                                idxs=np.concatenate(idxs))
            
            numberings, idxs, seq_counts = [], [], 0

    if seq_counts > 0:
                                
        unique_id = str(uuid.uuid4())
        save_file = os.path.join(save_folder, "data-subset-{}-{}.npz".format(suffix, unique_id))                        
        np.savez_compressed(save_file, 
                            numberings=np.concatenate(numberings), 
                            idxs=np.concatenate(idxs))
        
    [os.remove(fname) for fname in list_of_files]