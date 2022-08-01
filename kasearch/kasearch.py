import os
import glob
import time
import subprocess
import requests

import numpy as np
from concurrent.futures.thread import ThreadPoolExecutor

from kasearch.identity_calculations import get_n_most_identical, slow_get_n_most_identical
from kasearch.meta_extract import ExtractMetadata

class SearchDB(ExtractMetadata):
    def __init__(self, 
                 database_path='oasdb-small', 
                 allowed_chain='Any', 
                 allowed_species='Any', 
                 n_jobs=None,
                 id_to_study_file=None,
                ):
        super().__init__()
        if database_path == 'oasdb-small': # Download a small version of OAS
            db_folder = os.path.join(os.path.dirname(__file__), "oasdb")
            os.makedirs(db_folder, exist_ok = True)
            
            if not glob.glob(os.path.join(db_folder, "oasdb_small*")):
                print("Downloading a small version of OAS (4.4GB) ...")
                
                url = "https://zenodo.org/record/6668747/files/oasdb_small.tar"
                tmp_file = os.path.join(db_folder, "tmp.tar")

                with open(tmp_file,'wb') as f: f.write(requests.get(url).content)
                
                subprocess.run(["tar", "-xf", tmp_file, "-C", db_folder], check = True) 
                os.remove(tmp_file)
               
            database_path = glob.glob(os.path.join(db_folder, "oasdb_small*"))[0]
        
        self.database_path = database_path
        self.n_jobs = n_jobs
        
        self._set_id_to_study(os.path.join(database_path, "id_to_study.txt"))
        self.__set_allowed_files(allowed_chain, allowed_species)
        self.__reset_current_best()
    
    def __reset_current_best(self):
        
        self._current_target_numbering = None
        self._current_target_ids = None
        self.current_best_identities = np.zeros((1, 3), np.float16) - 1
        self.current_best_ids = np.zeros((1, 3, 2), np.int32) - 1

    def __set_allowed_files(self, allowed_chain, allowed_species):
        
        self.allowed_normal_files = []
        self.allowed_unusual_files = []
        
        if allowed_species == 'Any': allowed_species = '*'
        allowed_species = [allowed_species] if isinstance(allowed_species, str) else allowed_species
        if allowed_chain == 'Any': allowed_chain = '*'
        
        for species in allowed_species:
            self.allowed_normal_files += glob.glob(os.path.join(self.database_path, 
                                                         allowed_chain, 
                                                         species, 
                                                         "*data-subset-normal-*.npz"))
            self.allowed_unusual_files += glob.glob(os.path.join(self.database_path, 
                                                               allowed_chain, 
                                                               species, 
                                                               "*data-subset-unusual-*.npz"))

    def __update_best(self, query, keep_best_n):
        chunk_best_identities, chunk_best_ids = get_n_most_identical(query,
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

        data_loader = DataLoader(self.allowed_normal_files[0])
        self._current_target_numbering = data_loader.data['numberings']
        self._current_target_ids = data_loader.data['idxs']
        
        for file in self.allowed_normal_files[1:]:
            data_loader = DataLoader(file)
            self.__update_best(query, keep_best_n)
            self._current_target_numbering = data_loader.data['numberings']
            self._current_target_ids = data_loader.data['idxs']
                
        self.__update_best(query, keep_best_n)
        
    def get_meta(self, n_query = 0, n_region = 0, n_sequences = 'all', n_jobs=1):
        
        if n_sequences == 'all':
            n_sequences = len(self.current_best_ids)
            
        assert n_query >= 0
        assert n_region >= 0
        assert n_sequences > 0
        
        metadf = self._extract_meta(self.current_best_ids[:n_sequences, n_region], n_jobs=n_jobs)
        metadf['Identity'] = self.current_best_identities[:n_sequences, n_region]
        return metadf
        
        
        
class DataLoader():
    def __init__(self, file):
        _pool = ThreadPoolExecutor()
        self.__data = _pool.submit(self.__load_data, file)
        
    def __load_data(self, file):
        data = np.load(file, allow_pickle=True)
        output = {}
        output['numberings'] = data['numberings']
        output['idxs'] = data['idxs']
        
        return output
        
    @property
    def data(self):
        while self.__data.running():  # Want to acess only if it has finished
            time.sleep(1)
        return self.__data.result()
    
