import os
import uuid
from dataclasses import dataclass
from multiprocessing import Pool

import numpy as np

from kasearch.species_anarci import number
from kasearch.canonical_alignment import canonical_alignment, canonical_alignment_oas, canonical_numbering_len

class PrepareDB:
    
    def __init__(self, db_path, oas_source=False):
        
        os.makedirs(db_path, exist_ok=True)
        self.db_path = db_path
        self._file_suffix = 0
        self._oas_source = oas_source
        
        self._abnormal_sequence = np.zeros(canonical_numbering_len, np.int8)
        self.tmpDB = tmpDB({}, {}, {}, {}, {})
    
    def _canonical_alignment_oas(self, sequence):
        
        try:
            return canonical_alignment_oas(sequence)
        except Exception:
            return self._abnormal_sequence
    
    def _canonical_alignment(self, sequence):
        try:
            return canonical_alignment(sequence)
        except Exception:
            return self._abnormal_sequence
    
    def _many_canonical_alignment(self, sequence_dataset):
        
        if self._oas_source:
            return np.array([self._canonical_alignment_oas(sequence) for sequence in sequence_dataset])
        else:
            return np.array([self._canonical_alignment(sequence) for sequence in sequence_dataset])
        
    def _update_tmpDB(self, sequences, sequence_alignments, sequence_idxs, chain, species):
        
        if chain not in self.tmpDB.sequences_count:
            self.tmpDB.sequences_count[chain] = {}
            self.tmpDB.sequences[chain] = {}
            self.tmpDB.sequences_idxs[chain] = {}
            self.tmpDB.abnormal_sequences[chain] = {}
            self.tmpDB.abnormal_sequences_idxs[chain] = {}
            
        if species not in self.tmpDB.sequences_count[chain]:
            self.tmpDB.sequences_count[chain][species] = 0
            self.tmpDB.sequences[chain][species] = []
            self.tmpDB.sequences_idxs[chain][species] = []
            self.tmpDB.abnormal_sequences[chain][species] = []
            self.tmpDB.abnormal_sequences_idxs[chain][species] = []
            
        self.tmpDB.sequences_count[chain][species] += len(sequences)
        
        abnormal_sequence_idxs = np.where(0 == sequence_alignments.sum(-1))
        normal_sequence_idxs = np.where(0 != sequence_alignments.sum(-1))

        if normal_sequence_idxs[0] != []:
            self.tmpDB.sequences[chain][species].append(sequence_alignments[normal_sequence_idxs])
            self.tmpDB.sequences_idxs[chain][species].append(sequence_idxs[normal_sequence_idxs])
        
        if abnormal_sequence_idxs[0] != []:
            self.tmpDB.abnormal_sequences[chain][species].append(
                sequences[abnormal_sequence_idxs])
            self.tmpDB.abnormal_sequences_idxs[chain][species].append(
                sequence_idxs[abnormal_sequence_idxs])
    
    def _save_data_subset(self, sequence_alignments, sequence_idxs, chain, species, suffix = ''):
        
        if suffix == None:
            suffix = self._file_suffix
            
        unique_id = str(uuid.uuid4())
        save_folder = os.path.join(self.db_path, chain, species)
        save_file = os.path.join(save_folder, "data-subset-{}-{}.npz".format(suffix, unique_id))
        os.makedirs(save_folder, exist_ok=True)
        
        np.savez_compressed(save_file, 
                            numberings=np.concatenate(sequence_alignments[chain][species]), 
                            idxs=np.concatenate(sequence_idxs[chain][species]))
        
    def prepare_database(self, sequence_dataset, file_id, sequence_lines = None, chain='Heavy', species='Human'):
        """
        Prepares the new database. If a subset contains more than 50 million sequences, save that subset.
        """
        
        if not sequence_lines:
            sequence_lines = range(len(sequence_dataset))
        
        sequence_alignments = self._many_canonical_alignment(sequence_dataset)
        sequence_idxs = np.array([[file_id, i] for i in sequence_lines], np.int32)
        
        self._update_tmpDB(sequence_dataset, sequence_alignments, sequence_idxs, chain, species)
        
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
                    
                if self.tmpDB.abnormal_sequences[chain][species] != []:
                    self._save_data_subset(self.tmpDB.abnormal_sequences, 
                                           self.tmpDB.abnormal_sequences_idxs, 
                                           chain, species,
                                           suffix = 'unusual'
                                          )
    
    
    
def flatten(xss):
    return [x for xs in xss for x in xs]

def unflatten(xss):
    return [xss[x:x+2] for x in range(0, len(xss), 2)]
    
@dataclass   
class tmpDB: 
    sequences_count: dict
    sequences: dict
    sequences_idxs: dict
    abnormal_sequences: dict
    abnormal_sequences_idxs: dict