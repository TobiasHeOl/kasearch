from dataclasses import dataclass
from multiprocessing import Pool

from kasearch.species_anarci import number
from kasearch.canonical_alignment import canonical_alignment

class SetQueries:
    """
    Prepares queries for KA-Search. 
    """
    def __init__(self, queries, allowed_species=None, n_jobs=None):
        
        self.n_jobs = n_jobs
        self.allowed_species = [i.lower() for i in allowed_species]
        self.queries = self.__set_queries(queries)
        
        
    def _set_query(self, query):
        
        sequence = query
        
        numbered_query, chain = number(sequence, allowed_species=self.allowed_species)
        aligned_query = canonical_alignment(numbered_query)
        
        return Query(aligned_query, chain)
    
    def __set_queries(self, queries):
        size = len(queries)

        with Pool(processes=self.n_jobs) as pool:
            return pool.map(self._set_query, queries, chunksize=size // self.n_jobs)
        
        
@dataclass   
class Query:
    
    aligned_query: None
    chain: str
