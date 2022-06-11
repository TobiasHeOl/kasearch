import os
import ast
import pkg_resources

import pandas as pd


class ExtractMetadata:
    def __init__(self):
        
        with open(pkg_resources.resource_filename(__name__, "specific_studies.txt"), "r") as handle:
            self.specific_studies = ast.literal_eval(handle.readlines()[0])
            
        with open(pkg_resources.resource_filename(__name__, "public_studies.txt"), "r") as handle:
            self.public_studies = ast.literal_eval(handle.readlines()[0])
            
    def get_meta(self, study_id, line_id):
        
        study_file = self.public_studies[self.specific_studies[study_id]]
        return pd.read_csv(study_file, header=1, skiprows=list(range(2,line_id+2)), nrows=1)