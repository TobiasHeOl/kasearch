from multiprocessing import Pool

import numpy as np

from kasearch.species_anarci import number
from kasearch.canonical_alignment import canonical_alignment, canonical_alignment_oas, canonical_numbering_len

class AlignSequences:
    """
    Align sequences for KA-Search. 
    """
    def __init__(self, allowed_species=['Human', 'Mouse'], n_jobs=1, oas_source=False):
        
        self._oas_source = oas_source
        self.n_jobs = n_jobs
        self._abnormal_sequence = np.zeros(canonical_numbering_len, np.int8)
        
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
        
        else:
            n_jobs = len(seqs) if len(seqs) <  self.n_jobs else self.n_jobs
            chunksize=len(seqs) // n_jobs

            with Pool(processes=n_jobs) as pool:
                return np.array(pool.map(self._canonical_alignment, seqs, chunksize=chunksize))
        
    def __call__(self, seqs):
        
        if isinstance(seqs, str): seqs = [seqs]
        
        return self._many_canonical_alignment(seqs)

