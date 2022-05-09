import numpy as np
import os
import glob
from multiprocessing import Process
from kasearch.identity_calculations import get_n_most_identical, slow_get_n_most_identical


class SearchOAS:
    def __init__(self, database_path, chain='Heavy', allowed_species=None, n_jobs=None):

        self.database_path = database_path
        self.chain = chain
        self.n_jobs = n_jobs
        
        self.__set_allowed_files(allowed_species)
        self.current_target_numbering = None
        self.current_target_ids = None
        self.__reset_current_best()
    
    def __reset_current_best(self):
        
        self.current_best_identities = np.zeros((1, 3), np.float16) - 1
        self.current_best_ids = np.zeros((1, 3, 2), np.int32) - 1

    def __set_allowed_files(self, allowed_species):
        self.allowed_files = []
        for species in allowed_species:
            self.allowed_files += glob.glob(os.path.join(self.database_path, self.chain, species, "*.npz"))
        
    def __load_data_chunk(self, path=None):
        if path is not None:
            data = np.load(path)
            self.current_target_numbering = data['numberings']
            self.current_target_ids = data['idxs']

    def __update_best(self, query, keep_best_n):
        
        chunk_best_identities, chunk_best_ids = get_n_most_identical(query.aligned_query,
                                                                     self.current_target_numbering,
                                                                     self.current_target_ids, 
                                                                     n=keep_best_n,
                                                                     n_jobs=self.n_jobs)

        all_identities = np.concatenate([chunk_best_identities, self.current_best_identities])
        all_ids = np.concatenate([chunk_best_ids, self.current_best_ids])

        order = np.argsort(-all_identities, axis=0)
        self.current_best_identities = np.take_along_axis(all_identities, order, axis=0)[:keep_best_n]
        self.current_best_ids = np.take_along_axis(all_ids, order[:, :, None], axis=0)[:keep_best_n]
        
    def search(self, query, keep_best_n=10, reset_best=True):
        
        if reset_best == True:
            self.__reset_current_best()
        
        self.allowed_files = self.allowed_files[:1]
        
        for file in self.allowed_files:
            self.__load_data_chunk(file)
            self.__update_best(query, keep_best_n)
