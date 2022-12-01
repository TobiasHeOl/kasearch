import os
import glob
import subprocess
import requests

class InitiateDatabase:
    def __init__(self):
        pass
    
    def _set_database_path(self, database_path = 'oasdb-tiny'):
        """
        Sets database path. If none given, downloads a small version of OAS (4.4GB).
        """
        
        if database_path == 'oasdb-small': # Download a small version of OAS
            database_path = download_small_oas()
            
        if database_path == 'oasdb-tiny': # Download a tiny version of OAS
            database_path = download_tiny_oas()
        
        self.database_path = database_path
        
    def _set_files_to_search(self, allowed_chain, allowed_species):
        """
        Sets files to search.
        """
        
        self.files_to_search_normal = []
        self.files_to_search_unusual = []
        
        if allowed_species == 'Any': allowed_species = '*'
        allowed_species = [allowed_species] if isinstance(allowed_species, str) else allowed_species
        if allowed_chain == 'Any': allowed_chain = '*'
        
        for species in allowed_species:
            self.files_to_search_normal += glob.glob(os.path.join(self.database_path, 
                                                         allowed_chain, 
                                                         species, 
                                                         "*data-subset-normal-*.npz"))
            self.files_to_search_unusual += glob.glob(os.path.join(self.database_path, 
                                                               allowed_chain, 
                                                               species, 
                                                               "*data-subset-unusual-*.npz"))
            
            
def download_small_oas():
    
    db_folder = os.path.join(os.path.dirname(__file__), "oasdb")
    os.makedirs(db_folder, exist_ok = True)

    if not glob.glob(os.path.join(db_folder, "oasdb_small*")):
        print("Downloading a small version of OAS (4.4GB) ...")
        url = "https://zenodo.org/record/7384311/files/oasdb_small.tar"
        tmp_file = os.path.join(db_folder, "tmp.tar")

        with open(tmp_file,'wb') as f: f.write(requests.get(url).content)

        subprocess.run(["tar", "-xf", tmp_file, "-C", db_folder], check = True) 
        os.remove(tmp_file)

    return glob.glob(os.path.join(db_folder, "oasdb_small*"))[0]
    
def download_tiny_oas():
    
    db_folder = os.path.join(os.path.dirname(__file__), "oasdb")
    os.makedirs(db_folder, exist_ok = True)

    if not glob.glob(os.path.join(db_folder, "oas-aligned-tiny*")):
        print("Downloading a tiny version of OAS (~400MB) ...")

        url = "https://zenodo.org/record/7384311/files/oas-aligned-tiny.tar"
        tmp_file = os.path.join(db_folder, "tmp.tar")

        with open(tmp_file,'wb') as f: f.write(requests.get(url).content)

        subprocess.run(["tar", "-xf", tmp_file, "-C", db_folder], check = True) 
        os.remove(tmp_file)

    return glob.glob(os.path.join(db_folder, "oas-aligned-tiny*"))[0]