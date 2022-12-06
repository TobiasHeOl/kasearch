from multiprocessing import Pool

import numpy as np

from kasearch.species_anarci import number, many_number
from kasearch.canonical_alignment import canonical_alignment, canonical_alignment_oas, canonical_numbering_len

class AlignSequences:
    """
    Canonical alignment of sequences for KA-Search.  
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
        """Canonical alignment of a single sequence.
        
        Parameters
        ----------
        seq : str
            Antibody sequence to align

        Returns
        -------
        numpy array
            Canonical alignment of a single sequence
        """ 
        
        try:
            numbered_sequence, _ = number(seq, allowed_species=self.allowed_species)
        except Exception:
            raise ValueError(f"Sequence can not be numbered by ANARCI. This may be not an antibody variable domain.\nSequence that breaks: {seq}")
        try:
            return canonical_alignment(numbered_sequence)
        except Exception:
            raise ValueError(f"Sequence can not be aligned with the canonical alignment.\nSequence that breaks: {seq}")
        
    def _canonical_alignment_oas(self, seq):
        """Canonical alignment of OAS derived ANARCI numberings.
        
        Parameters
        ----------
        seq : dict
            OAS derived ANARCI numberings

        Returns
        -------
        numpy array
            Canonical alignment of a single sequence
        """
        
        try:
            return canonical_alignment_oas(seq)
        except Exception:
            return self._unusual_sequence
        
    def _fast_canonical_alignment(self, numbered_seq):
        """Canonical alignment of ANARCI numberings.
        
        Parameters
        ----------
        seq : dict
            Nested ANARCI numberings

        Returns
        -------
        numpy array
            Canonical alignment of a single sequence
        """
        
        try:
            return canonical_alignment(numbered_seq[0][0])
        except Exception:
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
        
        if self._oas_source:
            return np.array([self._canonical_alignment_oas(seq) for seq in seqs])
        
        elif self._if_fast:
            numbered_seqs = many_number(seqs, allowed_species=self.allowed_species, n_jobs=self.n_jobs)
            assert numbered_seqs, "Target DB contains sequences breaking ANARCI."
            
            n_jobs = len(numbered_seqs) if len(numbered_seqs) <  self.n_jobs else self.n_jobs
            chunksize=len(numbered_seqs) // n_jobs

            with Pool(processes=n_jobs) as pool:
                return np.array(pool.map(self._fast_canonical_alignment, numbered_seqs, chunksize=chunksize))
        
        else:
            n_jobs = len(seqs) if len(seqs) <  self.n_jobs else self.n_jobs
            chunksize=len(seqs) // n_jobs

            with Pool(processes=n_jobs) as pool:
                return np.array(pool.map(self._canonical_alignment, seqs, chunksize=chunksize))
        
    def __call__(self, seqs):
        
        if isinstance(seqs, str): seqs = [seqs]
        
        return self._many_canonical_alignment(seqs)

