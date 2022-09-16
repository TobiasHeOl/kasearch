from multiprocessing import Pool

import numpy as np

from kasearch.species_anarci import number, many_number
from kasearch.canonical_alignment import canonical_alignment, canonical_alignment_oas, canonical_numbering_len

class AlignSequences:
    """
    Align sequences for KA-Search. 
    """
    def __init__(self, allowed_species=['Human', 'Mouse'], n_jobs=1, oas_source=False, if_fast=False):
        
        self._oas_source = oas_source
        self._if_fast = if_fast
        self.n_jobs = n_jobs
        self._unusual_sequence = np.zeros(canonical_numbering_len, np.int8)
        
        if allowed_species:
            self.allowed_species = [i.lower() for i in allowed_species]
        else:
            self.allowed_species = None
        
    def _canonical_alignment(self, seq):
        """
        Alignment of a single sequence.
        """       
        try:
            numbered_sequence, _ = number(seq, allowed_species=self.allowed_species)
            return canonical_alignment(numbered_sequence)
        except Exception:
            return self._unusual_sequence
        
    def _canonical_alignment_oas(self, seq):
        """
        Alignment of ANARCI results from OAS.
        """
        try:
            return canonical_alignment_oas(seq)
        except Exception:
            return self._unusual_sequence
            
    def _many_canonical_alignment(self, seqs):
        """
        Alignment of multiple sequences.
        """
        if self._oas_source:
            return np.array([self._canonical_alignment_oas(seq) for seq in seqs])
        
        elif self._if_fast:
            numbered_seqs = many_number(seqs, allowed_species=self.allowed_species, n_jobs=self.n_jobs)
            assert numbered_seqs, "Target DB contains sequences breaking ANARCI."
                        
            return np.array([canonical_alignment(seq[0][0]) if seq!=None else self._unusual_sequence for seq in numbered_seqs])
        
        else:
            n_jobs = len(seqs) if len(seqs) <  self.n_jobs else self.n_jobs
            chunksize=len(seqs) // n_jobs

            with Pool(processes=n_jobs) as pool:
                return np.array(pool.map(self._canonical_alignment, seqs, chunksize=chunksize))
        
    def __call__(self, seqs):
        
        if isinstance(seqs, str): seqs = [seqs]
        
        return self._many_canonical_alignment(seqs)

