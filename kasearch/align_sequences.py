from multiprocessing import Pool

import numpy as np

from kasearch.anarci_numbering import number_many_at_once
from kasearch.canonical_alignment import canonical_alignment, canonical_alignment_oas, canonical_numbering_len

class AlignSequences:
    """
    Canonical alignment of sequences for KA-Search.  
    """
    
    def __init__(self, allowed_species=['Human', 'Mouse'], n_jobs=1, from_oas=False, strict=True, fast_implementation=False):
        
        self._from_oas = from_oas
        self._fast_implementation= fast_implementation
        self.strict = strict
        self.n_jobs = n_jobs
        self._unusual_sequence = np.zeros(canonical_numbering_len, np.int8)
        
        if allowed_species:
            self.allowed_species = [i.lower() for i in allowed_species]
        else:
            self.allowed_species = None
        
    def _canonical_alignment(self, numbered_seq):
        """Canonical alignment of ANARCI numberings.
        
        Parameters
        ----------
        seq : dict
            Either OAS derived ANARCI numberings or nested ANARCI numberings

        Returns
        -------
        numpy array
            Canonical alignment of a single sequence
        """
        
        try:
            if self._from_oas:
                return canonical_alignment_oas(numbered_seq)
            else:
                return canonical_alignment(numbered_seq)
        except Exception:
            
            if self.strict:
                raise ValueError(f"At least one sequence cannot be aligned with the canonical alignment.")
            else:
                return self._unusual_sequence
            
    def _many_canonical_alignment(self, seqs):
        """Canonical alignment of many sequences.
        
        Parameters
        ----------
        seqs : list
            List of antibody sequences

        Returns
        -------
        numpy array
            Canonical alignments of many sequences
        """
        
        if self._from_oas:
            numbered_seqs = seqs
        else:
            numbered_seqs = number_many_at_once(seqs, allowed_species=self.allowed_species, strict=self.strict, ncpu=self.n_jobs)
        
        
        if self.n_jobs == 1:
            return np.array([self._canonical_alignment(numbered_seq) for numbered_seq in numbered_seqs])
        else:
            n_jobs = min(len(numbered_seqs), self.n_jobs)
            chunksize = len(numbered_seqs) // n_jobs

            with Pool(processes=n_jobs) as pool:
                return np.array(pool.map(self._canonical_alignment, numbered_seqs, chunksize=chunksize))
        
    def __call__(self, seqs):
        
        if isinstance(seqs, str): seqs = [seqs]
        
        return self._many_canonical_alignment(seqs)

