import uuid
import glob
import os

import numpy as np


def merge_files(data_folder, data_file_size = 5_000_000):
    """
    Merges the files into files containing "data_file_size" of sequences. Default is 5 million.
    """
    for subfolder in glob.glob(os.path.join(data_folder, '*', '*')):

        normal_files = glob.glob(os.path.join(subfolder, 'data-subset-normal-*.npz'))
        unusual_files = glob.glob(os.path.join(subfolder, 'data-subset-unusual-*.npz'))

        merge_subfolder(list_of_files=normal_files, save_folder=subfolder,  data_file_size=data_file_size, suffix='normal')
        merge_subfolder(list_of_files=unusual_files, save_folder=subfolder, data_file_size=data_file_size, suffix='unusual')

def chunks(lst, n):
    """
    Yield successive n-sized chunks from lst.
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
        
def merge_subfolder(list_of_files, save_folder, data_file_size = 5_000_000, suffix='-0-'): 
    """
    Merges the files in a subfolder into files containing "data_file_size" of sequences. Default is 50 million.
    """
    numberings, idxs, seq_counts = [], [], 0
    
    
    
    for data, file_name in [(np.load(fname, allow_pickle=True), fname) for fname in list_of_files]:
        numberings.append(data['numberings'])
        idxs.append(data['idxs'])
        seq_counts += data['idxs'].shape[0]
        
        if seq_counts > data_file_size:            
            for sub_numberings, sub_idxs in zip(chunks(np.concatenate(numberings), data_file_size), chunks(np.concatenate(idxs), data_file_size)):
                
                if sub_idxs.shape[0] == data_file_size:
                    
                    unique_id = str(uuid.uuid4())
                    save_file = os.path.join(save_folder, "data-subset-{}-{}.npz".format(suffix, unique_id))
                    np.savez_compressed(save_file, 
                                        numberings=sub_numberings, 
                                        idxs=sub_idxs)

            numberings, idxs, seq_counts = [], [], 0
            numberings.append(sub_numberings)
            idxs.append(sub_idxs)
            seq_counts += sub_idxs.shape[0]
        
        
        os.remove(file_name)
        
    if seq_counts > 0:
                                
        unique_id = str(uuid.uuid4())
        save_file = os.path.join(save_folder, "data-subset-{}-{}.npz".format(suffix, unique_id))                        
        np.savez_compressed(save_file, 
                            numberings=np.concatenate(numberings), 
                            idxs=np.concatenate(idxs))



    
    
        
        