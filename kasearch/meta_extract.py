import os
import ast
import json
import pkg_resources
from joblib import Parallel, delayed

import pandas as pd


class ExtractMetadata:
    def __init__(self, id_to_study_file=None):
        
        if id_to_study_file == None:
            id_to_study_file = pkg_resources.resource_filename(__name__, "id_to_study.txt")
        
        with open(id_to_study_file, "r") as handle:
            self.id_to_study = ast.literal_eval(handle.readlines()[0])
            
    def _get_single_meta(self, idx):

        study_id, line_id = idx
        study_file = self.id_to_study[study_id]
        sequence_meta = pd.Series(json.loads(','.join(pd.read_csv(study_file, nrows=0).columns)))
        sequence_data = pd.read_csv(study_file, 
                                    header=1, 
                                    skiprows=list(range(2,line_id+2)), 
                                    nrows=1).iloc[0]

        return pd.concat([sequence_data, sequence_meta])
    
    def get_meta(self, idxs, n_jobs=1):
        
        idxs = [idxs] if any(isinstance(i, int) for i in idxs) else idxs
        n_jobs = len(idxs) if len(idxs) <  n_jobs else n_jobs
        chunksize=len(idxs) // n_jobs

        return pd.concat(Parallel(n_jobs=n_jobs)(delayed(self._get_single_meta)(idx) for idx in idxs), axis=1).T
 