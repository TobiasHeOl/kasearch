import os
import shutil
import uuid
import collections
from dataclasses import dataclass
from multiprocessing import Pool

import numpy as np
import pandas as pd

from kasearch import AlignSequences
from kasearch.merge_db import merge_files
    
    
class TemporaryDataHolder:
    
    def __init__(self):
        self.sequences_count = collections.defaultdict(dict)
        self.sequences = collections.defaultdict(dict)
        self.sequences_idxs = collections.defaultdict(dict)
        self.unusual_sequences = collections.defaultdict(dict)
        self.unusual_sequences_idxs = collections.defaultdict(dict)  
    
    def set_empty_nested_dict(self, chain, species):
        
        self.sequences_count[chain][species] = 0
        self.sequences[chain][species] = []
        self.sequences_idxs[chain][species] = []
        self.unusual_sequences[chain][species] = []
        self.unusual_sequences_idxs[chain][species] = []
        
    def add_to_prepared_data(self, chain, species, sequences, sequence_alignments, sequence_idxs):
    
        if species not in self.sequences_count[chain]:
            self.set_empty_nested_dict(chain, species)
            
        self.sequences_count[chain][species] += len(sequences)
        
        unusual_sequence_idxs = np.where(0 == sequence_alignments.sum(-1))
        normal_sequence_idxs = np.where(0 != sequence_alignments.sum(-1))

        if not isinstance(normal_sequence_idxs[0], list):
            self.sequences[chain][species].append(sequence_alignments[normal_sequence_idxs])
            self.sequences_idxs[chain][species].append(sequence_idxs[normal_sequence_idxs])
        
        if not isinstance(unusual_sequence_idxs[0], list):
            self.unusual_sequences[chain][species].append(
                sequences[unusual_sequence_idxs])
            self.unusual_sequences_idxs[chain][species].append(
                sequence_idxs[unusual_sequence_idxs])
            
    def save_data_subset(self, chain, species):
        
        save_folder = os.path.join(self.db_path, chain, species)
        os.makedirs(save_folder, exist_ok=True)    
        
        if self.sequences[chain][species] != []:
            np.savez_compressed(file = os.path.join(save_folder, f"data-subset-normal-{uuid.uuid4()}.npz"), 
                                numberings = np.concatenate(self.sequences[chain][species]), 
                                idxs = np.concatenate(self.sequences_idxs[chain][species]))
        
        if self.unusual_sequences[chain][species] != []:
            np.savez_compressed(file = os.path.join(save_folder, f"data-subset-unusual-{uuid.uuid4()}.npz"), 
                                numberings = np.concatenate(self.unusual_sequences[chain][species]), 
                                idxs = np.concatenate(self.unusual_sequences_idxs[chain][species]))
        
        self.set_empty_nested_dict(chain, species)
        
    def save_data_all(self):
        
        for chain, species_dict in self.sequences_idxs.items():
            for species in species_dict.keys():
                self.save_data_subset(chain, species)
                
            
class PrepareDB(TemporaryDataHolder):
    
    def __init__(self, db_path, n_jobs=1, from_scratch=False, from_oas=False):
        super().__init__()
        
        self.db_path, self.n_jobs, self.from_oas = db_path, n_jobs, from_oas
        
        if from_scratch and os.path.exists(db_path): shutil.rmtree(db_path)
        
        os.makedirs(db_path, exist_ok=True)
        self.id_to_study = {}
        
    def prepare_sequences(
        self, 
        sequence_file, 
        file_id, 
        chain = 'Heavy', 
        species ='Any', 
        seq_column_name = 'sequence_alignment_aa',
        strict = False,
        sequence_line_idx = None, 
        sequence_numberings = None,
    ):
        """
        Prepares a new database. 
        """
        
        os.makedirs(os.path.join(self.db_path, 'extra_data'), exist_ok=True)
        
        if self.from_oas:
            self.id_to_study[file_id] = sequence_file
        else:
            self.id_to_study[file_id] = os.path.basename(sequence_file)
        
        if sequence_numberings is None:
            sequences = pd.read_csv(sequence_file, header=1, usecols=[seq_column_name]).iloc[:,0].values
        else:
            sequences = sequence_numberings
        
        
        anarci_species = [species] if species != 'Any' else ['Human', 'Mouse']
            
        sequence_alignments = AlignSequences(
            n_jobs=self.n_jobs, 
            allowed_species=anarci_species,
            from_oas=self.from_oas, 
            strict=strict
        )(sequences)
        
        if not sequence_line_idx: sequence_line_idx = range(len(sequences))

        sequence_idxs = np.array([[file_id, i] for i in sequence_line_idx], np.int32)
        
        self.add_to_prepared_data(chain, species, sequences, sequence_alignments, sequence_idxs)
        
        if self.sequences_count[chain][species]>4_000_000:
            self.save_data_subset(chain, species)            
            
    def finalize_prepared_files(self, prepared_file_size = 5_000_000):
        
        self.save_data_all()
        
        merge_files(self.db_path, data_file_size = prepared_file_size)
        
        with open(os.path.join(self.db_path, "id_to_study.txt"), "w") as handle: 
            handle.write(str(self.id_to_study))
            
            
            
            
            
