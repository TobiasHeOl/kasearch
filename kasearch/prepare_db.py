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


class PrepareDB:
    
    def __init__(self, db_path, n_jobs=1, from_scratch=False, oas_source=False):
        
        if from_scratch: 
            if os.path.exists(db_path): shutil.rmtree(db_path)
        
        os.makedirs(db_path, exist_ok=True)
        self.db_path = db_path
        self.id_to_study = {}
        self.n_jobs = n_jobs
        self._oas_source = oas_source
        
        self.tmpDB = tmpDB(collections.defaultdict(dict), 
                           collections.defaultdict(dict), 
                           collections.defaultdict(dict), 
                           collections.defaultdict(dict), 
                           collections.defaultdict(dict))
        
    def _update_tmpDB(self, sequences, sequence_alignments, sequence_idxs, chain, species):
        
        if species not in self.tmpDB.sequences_count[chain]:
            self.tmpDB.sequences_count[chain][species] = 0
            self.tmpDB.sequences[chain][species] = []
            self.tmpDB.sequences_idxs[chain][species] = []
            self.tmpDB.unusual_sequences[chain][species] = []
            self.tmpDB.unusual_sequences_idxs[chain][species] = []
            
        self.tmpDB.sequences_count[chain][species] += len(sequences)
        
        unusual_sequence_idxs = np.where(0 == sequence_alignments.sum(-1))
        normal_sequence_idxs = np.where(0 != sequence_alignments.sum(-1))

        if normal_sequence_idxs[0] != []:
            self.tmpDB.sequences[chain][species].append(sequence_alignments[normal_sequence_idxs])
            self.tmpDB.sequences_idxs[chain][species].append(sequence_idxs[normal_sequence_idxs])
        
        if unusual_sequence_idxs[0] != []:
            self.tmpDB.unusual_sequences[chain][species].append(
                sequences[unusual_sequence_idxs])
            self.tmpDB.unusual_sequences_idxs[chain][species].append(
                sequence_idxs[unusual_sequence_idxs])
    
    def _save_data_subset(self, sequence_alignments, sequence_idxs, chain, species, suffix = ''):
 
        save_folder = os.path.join(self.db_path, chain, species)
        save_file = os.path.join(save_folder, "data-subset-{}-{}.npz".format(suffix, str(uuid.uuid4())))
        os.makedirs(save_folder, exist_ok=True)
        
        np.savez_compressed(save_file, 
                            numberings=np.concatenate(sequence_alignments[chain][species]), 
                            idxs=np.concatenate(sequence_idxs[chain][species]))
        
    def prepare_sequences(
        self, 
        sequence_file, 
        file_id, 
        chain='Heavy', 
        species='Human', 
        seq_column_name = 'sequence_alignment_aa',
        sequence_lines = None, 
        pre_calculated_anarci = None
    ):
        """
        Prepares a new database. 
        """
        
        os.makedirs(os.path.join(self.db_path, 'extra_data'), exist_ok=True)

        self.id_to_study[file_id] = os.path.basename(sequence_file)
        
        if pre_calculated_anarci is None:
            sequences = pd.read_csv(sequence_file, header=1, usecols=[seq_column_name]).iloc[:,0].values
        else:
            sequences = pre_calculated_anarci
        
        if species != 'Any':
            anarci_species = [species]
        else:
            anarci_species = ['Human', 'Mouse']
            
        sequence_alignments = AlignSequences(
            n_jobs=self.n_jobs, 
            allowed_species=anarci_species,
            oas_source=self._oas_source, 
            if_fast=True
        )(sequences)
        
        if not sequence_lines: sequence_lines = range(len(sequences))

        sequence_idxs = np.array([[file_id, i] for i in sequence_lines], np.int32)
        
        self._update_tmpDB(sequences, sequence_alignments, sequence_idxs, chain, species)
        
        if self.tmpDB.sequences_count[chain][species]>4_000_000:
            
            self._save_data_subset(self.tmpDB.sequences, 
                                   self.tmpDB.sequences_idxs, 
                                   chain, species,
                                   suffix = 'normal'
                                  )
            
            self.tmpDB.sequences_count[chain][species] = 0
            self.tmpDB.sequences[chain][species] = []
            self.tmpDB.sequences_idxs[chain][species] = []
            
    
       
    def save_database(self):
        
        for chain, species_dict in self.tmpDB.sequences_idxs.items():
            for species in species_dict.keys():
                if self.tmpDB.sequences[chain][species] != []:
                    self._save_data_subset(self.tmpDB.sequences, 
                                           self.tmpDB.sequences_idxs, 
                                           chain, species,
                                           suffix = 'normal'
                                          )
                    
                if self.tmpDB.unusual_sequences[chain][species] != []:
                    self._save_data_subset(self.tmpDB.unusual_sequences, 
                                           self.tmpDB.unusual_sequences_idxs, 
                                           chain, species,
                                           suffix = 'unusual'
                                          )
                    
        with open(os.path.join(self.db_path, "id_to_study.txt"), "w") as handle: 
            handle.write(str(self.id_to_study))
            
    def merge_sequence_files(self):
        
        merge_files(self.db_path)
    
    
    
def flatten(xss):
    return [x for xs in xss for x in xs]

def unflatten(xss):
    return [xss[x:x+2] for x in range(0, len(xss), 2)]
    
@dataclass   
class tmpDB: 
    sequences_count: dict
    sequences: dict
    sequences_idxs: dict
    unusual_sequences: dict
    unusual_sequences_idxs: dict