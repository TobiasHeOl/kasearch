import os, glob, json, random

from multiprocessing import Pool

import pandas as pd
import numpy as np
import uuid
import shutil

from kasearch import PrepareDB
from kasearch.merge_db import merge_files

    
class PrepareOASdb:
    """
    Function used to prepare OAS data for KA-Search. It extracts ANARCI numberings from OAS files, 
    transforms the numberings into encodings, and finally merging the encodings into files of a certain size.
    
    Additionally, it also creates a separate file for unusual sequences, i.e. sequences which have positions 
    other than the 200 used for the fast calculations. 
    """
    
    def __init__(
        self, 
        final_db_folder, 
        local_oas_path,
        filter_data = False,
        n_jobs = 20,
        public_path = "http://opig.stats.ox.ac.uk/webapps/ngsdb/"
    ):
        self.public_path = public_path
        self.n_jobs = n_jobs
        self.final_db_folder = final_db_folder
        self.filter_data = filter_data
        if os.path.exists(final_db_folder): shutil.rmtree(final_db_folder)
        os.makedirs(final_db_folder, exist_ok=True)

        self.set_id_to_study(local_oas_path)
        
        
    def set_id_to_study(self, local_oas_path):
        """
        Creates the id_to_study dictionary. This dictionary is used for finding the relevant 
        meta data for a given id, using the get_meta function

        Additionally, it also adds an index to the data_unit_files, required for multiprocessing
        """
        data_unit_files = glob.glob(os.path.join(local_oas_path, 'unpaired','*/*/*.csv.gz'))
        print('# of unpaired files:', len(data_unit_files))

        self.data_unit_files = [[num, file] for num, file in enumerate(data_unit_files)]       
        
        id_to_study = {}
        for num, file in enumerate(data_unit_files):
            id_to_study[num] = os.path.join(self.public_path, file.strip(local_oas_path))
            
        self.id_to_study = id_to_study
        
        with open(os.path.join(self.final_db_folder, "id_to_study.txt"), "w") as handle:
            handle.write(str(self.id_to_study))
        

    def process_data(self, data_files, sequence_lines=None):
        """
        Wrapper around PrepareDB, for ANARCI numberings from OAS. Optimized for faster processing.
        
        Additionally, it also allows for filtering to create a smaller and cleaner version to search.
        """

        preparedb = PrepareDB(db_path=self.final_db_folder, n_jobs=1, from_oas=True, from_scratch=False)
        for num, data_file in unflatten(data_files):

            metadata = json.loads(','.join(pd.read_csv(data_file, nrows=0).columns))
            species, chain = format_species(metadata['Species']), metadata['Chain']

            DB = pd.read_csv(data_file, header=1, usecols=['ANARCI_numbering', 'ANARCI_status', 'Redundancy', 'sequence_alignment_aa'])

            if self.filter_data:
                DB.query("Redundancy>=5", inplace=True)
                DB.query("~ANARCI_status.str.contains('Cysteine')", engine='python', inplace=True)
                DB.query("~ANARCI_status.str.contains('Unusual residue')", engine='python', inplace=True) 
                DB.query("~sequence_alignment_aa.str.contains('\*')", engine='python', inplace=True)
                DB.query("~sequence_alignment_aa.str.contains('X')", engine='python', inplace=True) 
                sequence_lines = list(DB.index)

                if DB.empty: continue # Otherwise prepare_sequences breaks

            preparedb.prepare_sequences(
                data_file, 
                file_id=num, 
                sequence_line_idx = sequence_lines, 
                chain=chain, 
                species=species,
                sequence_numberings=DB.ANARCI_numbering.values
            )

        preparedb.save_data_all()

        
    def process_many_files(self):
        """
        Function for initiating multiprocessing
        """
        
        data_unit_files = self.data_unit_files
        
        random.shuffle(data_unit_files)

        subset_size = min(len(data_unit_files), 50)
        data_file_subsets = [flatten(data_unit_files[x:x+subset_size]) for x in range(0, len(data_unit_files), subset_size)]

        n_jobs = len(data_file_subsets) if len(data_file_subsets) <  self.n_jobs else self.n_jobs

        with Pool(processes=n_jobs) as pool:
            return pool.map(self.process_data, data_file_subsets)
        


    def __call__(self, data_file_size = 50_000_000):

        self.process_many_files()

        merge_files(self.final_db_folder, data_file_size = data_file_size) # Merge folders into sets of 50 million sequences
    
        

    
def format_species(species):
    """
    This is not strictly needed, but results in a nicer formatted dataset
    """
    
    if 'his' in species.lower():
        return 'Humanised'
    elif 'kymouse' in species.lower():
        return 'Humanised'
    elif 'human' in species.lower():
        return 'Human'
    elif 'mouse' in species.lower():
        return 'Mouse'
    elif 'rat' in species.lower():
        return 'Rat'
    elif 'rabbit' in species.lower():
        return 'Rabbit'
    elif 'camel' in species.lower():
        return 'Camel'
    elif 'rhesus' in species.lower():
        return 'Rhesus'
    else:
        return 'Unknown'                                      
                                      
def flatten(xss):
    return [x for xs in xss for x in xs]

                                      
def unflatten(xss):
    return [xss[x:x+2] for x in range(0, len(xss), 2)]

    
    
def prepare_tiny_oas(final_db_folder, oasdb_small_folder):
    """
    Codebase for creating tiny OAS, a clean set of full human antibody variable domains.
    """
    
    
    if os.path.exists(final_db_folder): shutil.rmtree(final_db_folder)
    os.makedirs(final_db_folder, exist_ok=True)
    os.makedirs(os.path.join(final_db_folder, "extra_data"), exist_ok=True)
    os.makedirs(os.path.join(final_db_folder, "Heavy", "Human"), exist_ok=True)
    shutil.copyfile(os.path.join(oasdb_small_folder, "id_to_study.txt"), os.path.join(final_db_folder, "id_to_study.txt"))
    
    numberings = []
    idxs = []
    for file in glob.glob(os.path.join(oasdb_small_folder, "Heavy/Human/*normal*")):
        data = np.load(file, allow_pickle=True)
        numberings.append(data['numberings'])
        idxs.append(data['idxs']) 


    numberings = np.concatenate(numberings)
    idxs = np.concatenate(idxs)
    
    tmp_df = pd.DataFrame(numberings, columns=[f"col_{i}" for i in range(200)], dtype=np.int8).query("col_0!=0")
    print("OAS small size", tmp_df.shape)
    tmp_df = tmp_df.query("col_0!=124")
    print("# full abs", tmp_df.shape)

    tmp_df = tmp_df.drop_duplicates()
    print("# unique abs", tmp_df.shape)
    
    for set_of_idxs in np.array_split(tmp_df.index, 10):
        np.savez_compressed(
            os.path.join(final_db_folder, "Heavy", "Human", f"data-subset-normal-{uuid.uuid4()}.npz"), 
            numberings=numberings[set_of_idxs], 
            idxs=idxs[set_of_idxs]
        )