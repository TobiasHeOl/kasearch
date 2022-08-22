import os
import time

import numpy as np
from concurrent.futures.thread import ThreadPoolExecutor

from kasearch.identity_calculations import get_n_most_identical_multiquery, slow_get_n_most_identical
from kasearch.meta_extract import ExtractMetadata
from kasearch.initiate_db import InitiateDatabase
from kasearch.canonical_alignment import get_region_mask

class SearchDB(InitiateDatabase, ExtractMetadata):
    def __init__(self, 
                 database_path='oasdb-small', 
                 allowed_chain='Any', 
                 allowed_species='Any',
                 regions=['whole', 'cdrs', 'cdr3'],
                 length_matched=[False,True,True],
                 id_to_study_file=None,
                ):
        super().__init__()
        
        self.region_masks = np.stack([get_region_mask(region) for region in regions])
        self.length_matched = np.array(length_matched, dtype = bool)
        assert self.region_masks.shape[0] == self.length_matched.shape[0], "List of user-defined regions ({}) and 'if length match' ({})\
 are of different lengths. Please define a 'if length match' for each defined region.".format(self.region_masks.shape[0], self.length_matched.shape[0])
        
        self._set_database_path(database_path)
        self._set_files_to_search(allowed_chain, allowed_species)
        
        self._set_id_to_study(os.path.join(self.database_path, "id_to_study.txt"))
        
        assert len(self.files_to_search_normal) != 0, "DB does not contain data of {} chains from the {} species.".format(allowed_chain, allowed_species)
    
    def __reset_current_best(self, qsize=1):
        """
        Resets currently most similar sequences.
        """
        
        self._current_target_numbering = None
        self._current_target_ids = None
        self.current_best_identities = np.zeros((qsize, 1, self.region_masks.shape[0]), np.float16) - 1
        self.current_best_ids = np.zeros((qsize, 1, self.region_masks.shape[0], 2), np.int32) - 1

    def __update_best(self, query, keep_best_n):
        """
        Update the current most similar sequences.
        """
        
        chunk_best_identities, chunk_best_ids = get_n_most_identical_multiquery(query,
                                                                                self._current_target_numbering,
                                                                                self._current_target_ids, 
                                                                                n=keep_best_n,
                                                                                region_masks=self.region_masks, 
                                                                                length_matched=self.length_matched
                                                                               )
        
        all_identities = np.concatenate([chunk_best_identities, self.current_best_identities], axis=1)
        all_ids = np.concatenate([chunk_best_ids, self.current_best_ids], axis=1)

        order = np.argsort(-all_identities, axis=1)

        self.current_best_identities = np.take_along_axis(all_identities, order, axis=1)[:, :keep_best_n]
        self.current_best_ids = np.take_along_axis(all_ids, order[:, :, :, None], axis=1)[:, :keep_best_n]
        
    def search(self, query, keep_best_n=10):
        """
        Search database for sequences most similar to the queries.
        """
        
        self.__reset_current_best(query.shape[0])

        data_loader = DataLoader(self.files_to_search_normal[0])
        self._current_target_numbering = data_loader.data['numberings']
        self._current_target_ids = data_loader.data['idxs']
        
        for file in self.files_to_search_normal[1:]:
            data_loader = DataLoader(file)
            self.__update_best(query, keep_best_n)
            self._current_target_numbering = data_loader.data['numberings']
            self._current_target_ids = data_loader.data['idxs']
                
        self.__update_best(query, keep_best_n)
        
    def get_meta(self, n_query = 0, n_region = 0, n_sequences = 'all', n_jobs=1):
        """
        Get meta data for the current most similar sequences for a specific query and region.
        """
        
        if n_sequences == 'all':
            n_sequences = self.current_best_identities.shape[1]
            
        assert n_query >= 0
        assert n_region >= 0
        assert n_sequences > 0
        
        metadf = self._extract_meta(self.current_best_ids[n_query, :n_sequences, n_region], n_jobs=n_jobs)
        metadf['Identity'] = self.current_best_identities[n_query, :n_sequences, n_region]
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
    
