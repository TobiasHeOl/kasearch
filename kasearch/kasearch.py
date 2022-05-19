import numpy as np
import os
import glob
import concurrent
import time
from kasearch.identity_calculations import get_n_most_identical, slow_get_n_most_identical

class SearchOAS:
    def __init__(self, 
                 database_path, 
                 allowed_chain='Heavy', 
                 allowed_species=None, 
                 n_jobs=None):
        
        self.database_path = database_path
        self.n_jobs = n_jobs
        
        self.allowed_files = []
        self.__set_allowed_files(allowed_chain, allowed_species)
        
        self.__reset_current_best()
    
    def __reset_current_best(self):
        
        self._current_target_numbering = None
        self._current_target_ids = None
        self.current_best_identities = np.zeros((1, 3), np.float16) - 1
        self.current_best_ids = np.zeros((1, 3, 2), np.int32) - 1

    def __set_allowed_files(self, allowed_chain, allowed_species):
        
        if allowed_species == None: allowed_species = '*'      
        if allowed_chain == None: allowed_chain = '*'
        
        for species in allowed_species:
            self.allowed_files += glob.glob(os.path.join(self.database_path, allowed_chain, species, "*.npz"))

    def __update_best(self, query, keep_best_n):
        chunk_best_identities, chunk_best_ids = get_n_most_identical(query.aligned_query,
                                                                     self._current_target_numbering,
                                                                     self._current_target_ids, 
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

        #self.allowed_files = self.allowed_files[:2]
        data_loader = DataLoader(self.allowed_files[0])
        self._current_target_numbering = data_loader.data['numberings']
        self._current_target_ids = data_loader.data['idxs']
        
        for file in self.allowed_files[1:]:
            data_loader = DataLoader(file)
            self.__update_best(query, keep_best_n)
            self._current_target_numbering = data_loader.data['numberings']
            self._current_target_ids = data_loader.data['idxs']
                
        self.__update_best(query, keep_best_n)
        
        
class DataLoader():
    def __init__(self, file):
        _pool = concurrent.futures.ThreadPoolExecutor()
        self.__data = _pool.submit(self.__load_data, file)
        
    def __load_data(self, file):
        data = np.load(file)
        output = {}
        output['numberings'] = data['numberings']
        output['idxs'] = data['idxs']
        
        return output
        
    @property
    def data(self):
        while self.__data.running():  # Want to acess only if it has finished
            time.sleep(1)
        return self.__data.result()
    
