from dataclasses import dataclass
from multiprocessing import Pool

from kasearch.species_anarci import number
from kasearch.canonical_alignment import canonical_alignment

class SetQueries:
    """
    Prepares queries for KA-Search. 
    """
    def __init__(self, queries, allowed_species=None, n_jobs=1):
        
        self.n_jobs = n_jobs
        
        if allowed_species:
            self.allowed_species = [i.lower() for i in allowed_species]
        else:
            self.allowed_species = None
        
        self.queries = self.__set_queries(queries)
        
        
    def _set_query(self, query):
        
        sequence = query
        
        numbered_query, chain = number(sequence, allowed_species=self.allowed_species)
        aligned_query = canonical_alignment(numbered_query)
        
        return Query(aligned_query, chain)
    
    def __set_queries(self, queries):
        
        n_jobs = len(queries) if len(queries) <  self.n_jobs else self.n_jobs

        with Pool(processes=n_jobs) as pool:
            return pool.map(self._set_query, queries, chunksize=len(queries) // n_jobs)
        
        
@dataclass   
class Query:
    
    aligned_query: None
    chain: str
