
---

<div align="center">    
 
# KA-Search: Rapid and exhaustive sequence identity search of known antibodies


</div>


Antibodies with similar amino acid sequences, especially in the complementary-determining regions (CDRs), often share certain properties. It is often powerful to compare the sequence of an antibody of interest against natural antibody repertoires, as finding similar antibodies in nature can indicate likely specificity or immunogenicity. However, as the number of available antibody repertoire sequences has exceeded a billion and is continuing to grow, repertoire mining for highly similar sequences has become increasingly computationally expensive. Existing approaches are limited by either being low-throughput, non-exhaustive, not antibody-specific, or only searching against entire chain sequences. Therefore, there is a need for a specialized tool, optimized for a rapid and exhaustive search of any antibody region against all known antibodies, to better utilize the full number of available repertoire sequences.

Here, we introduce Known Antibody Search (KA-Search), a tool that allows for rapid search of the 1.7 billion antibodies in the Observed Antibody Space (OAS) database by sequence identity across either the whole chain, the CDRs, or a user defined antibody region. KA-Search can be used to find the most similar sequences from OAS within 20 minutes using 5 CPUs. We demonstrate how KA-Search can be used to obtain new insights about an antibody of interest. KA-Search is freely available at https://github.com/oxpig/kasearch.

-----------

# Install KA-Search

KA-Search is freely available and can be installed with pip.

~~~.sh
    pip install kasearch
~~~

or directly from github.

~~~.sh
    pip install -U git+https://github.com/oxpig/kasearch
~~~


Additionally, you need a version of [ANARCI](https://github.com/oxpig/ANARCI) in the same environment.

----------

# Searching with KA-Search

**A Jupyter notebook** showcasing KA-Search can be found [here](https://github.com/oxpig/kasearch/tree/main/examples). 

## Align query sequence

Align query sequences using the AlignSequences class.

```{r, engine='python', count_lines}

from kasearch import AlignSequences, SearchDB, ExtractMetadata, PrepareDB

raw_queries = [ 'QVKLQESGAELARPGASVKLSCKASGYTFTNYWMQWVKQRPGQGLDWIGAIYPGDGNTRYTHKFKGKATLTADKSSSTAYMQLSSLASEDSGVYYCARGEGNYAWFAYWGQGTTVTVSS', 'EVQLQQSGTVLARPGASVKMSCEASGYTFTNYWMHWVKQRPGQGLEWIGAIYPGNSDTSYIQKFKGKAKLTAVTSTTSVYMELSSLTNEDSAVYYCTLYDGYYVFAYWGQGTLVTVSA',
]

query_db = AlignSequences(raw_queries, # Sequences as strings to align.
                          n_jobs=1     # Allocated number for jobs/threads for the search.
                         )
query_db.db.aligned_seqs[0]
```
The alignment should look like this.

```console
array([[81, 86, 75,  0, 76, 81, 69, 83, 71, 65,  0, 69, 76, 65, 82, 80,
        71, 65, 83, 86, 75, 76, 83, 67, 75, 65, 83, 71, 89, 84, 70,  0,
         0,  0,  0,  0,  0,  0,  0,  0, 84, 78, 89, 87, 77, 81,  0, 87,
        86, 75, 81,  0, 82,  0, 80,  0, 71,  0, 81,  0,  0, 71,  0, 76,
        68,  0, 87, 73, 71, 65, 73, 89, 80, 71,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0, 68, 71, 78, 84, 82, 89,  0,  0, 84,  0,  0,
        72,  0,  0, 75, 70,  0,  0, 75,  0,  0,  0, 71, 75, 65, 84, 76,
        84, 65,  0, 68,  0,  0,  0, 75,  0, 83,  0,  0, 83, 83,  0,  0,
         0,  0, 84,  0, 65, 89, 77, 81, 76, 83, 83, 76, 65, 83,  0, 69,
        68, 83, 71, 86, 89, 89, 67, 65, 82, 71, 69, 71, 78,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0, 89, 65, 87, 70, 65, 89, 87, 71,  0, 81,
        71, 84, 84, 86, 84, 86, 83, 83],
       [69, 86, 81,  0, 76, 81, 81, 83, 71, 84,  0, 86, 76, 65, 82, 80,
        71, 65, 83, 86, 75, 77, 83, 67, 69, 65, 83, 71, 89, 84, 70,  0,
         0,  0,  0,  0,  0,  0,  0,  0, 84, 78, 89, 87, 77, 72,  0, 87,
        86, 75, 81,  0, 82,  0, 80,  0, 71,  0, 81,  0,  0, 71,  0, 76,
        69,  0, 87, 73, 71, 65, 73, 89, 80, 71,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0, 78, 83, 68, 84, 83, 89,  0,  0, 73,  0,  0,
        81,  0,  0, 75, 70,  0,  0, 75,  0,  0,  0, 71, 75, 65, 75, 76,
        84, 65,  0, 86,  0,  0,  0, 84,  0, 83,  0,  0, 84, 84,  0,  0,
         0,  0, 83,  0, 86, 89, 77, 69, 76, 83, 83, 76, 84, 78,  0, 69,
        68, 83, 65, 86, 89, 89, 67, 84, 76, 89, 68, 71, 89,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0, 89, 86, 70, 65, 89, 87, 71,  0, 81,
        71, 84, 76, 86, 84, 86, 83, 65]], dtype=int8)
```

## Initiate target database

Next, you initiate the target database to search against. Default is OAS, however, a description for creating your own database to search against can be found further down.

```{r, engine='python', count_lines}
DB_PATH = "../data/db-example"

oasdb = SearchDB(database_path=DB_PATH,  # DB path. Default will be to download a prepared version of OAS.
                 allowed_chain='Any',    # Search against a specific chain. Default is any chain.
                 allowed_species='Any',  # Search against a specific species. Default is any species.
                 n_jobs=5                # Allocated number for jobs/threads for the search.
                )
``` 

## Search database with aligned query

The initiated database can now be searched using an aligned query.
```{r, engine='python', count_lines}
oasdb.search(query_db.db.aligned_seqs[0], # Input can only be a single aligned sequence at a time.
             keep_best_n=2,               # You can define how many most similar sequences to return
             reset_best=True              # In cases where you want to search with the same sequence against multiple databases, you might want to not reset_best. (True is default)
            )

search_oas.current_best_identities
```

Returning the following best identities

```console
array([[0.79510413, 0.78571429, 0.83333333],
       [0.79161837, 0.78571429, 0.78409091]])
```

## Extract meta data from the N best results

The meta data of the N best matches can then be extracted using the N best indexes.

To get the meta data of the N best sequences from a specific region, the correct index needs to be specified.
- 0: whole sequence
- 1: CDRs
- 2: CDR3

e.g. "oasdb.current_best_ids\[:,0\]\" returns meta data for the N best whole sequence identities.

NB: The column "sequence_alignment_aa" holds the antibody sequence.


```{r, engine='python', count_lines}
meta_db = ExtractMetadata()

n_best_sequences = meta_db.get_meta(oasdb.current_best_ids[:,0])
n_best_sequences
```

Returns

```console
array([[0.79510413, 0.78571429, 0.83333333],
       [0.79161837, 0.78571429, 0.78409091]])
```

## Create custom database

To create your own database you first need to create a csv file in the OAS format. For an example file, look at data/custom-data-example.csv. This file consists of a dictionary containing the metadata in the first line and then rows of the individual sequences afterwards. Only the Species and Chain is strictly needed in the metadata, and only the amino acids sequence of the antibodies is required for each antibody sequence.

### 1. Format your data into OAS files

```{r, engine='python', count_lines}
import json
import os
import pandas as pd

from kasearch.species_anarci import number
from kasearch.merge_db import merge_files

metadata = {"Species":"Human", "Chain":"Heavy"}
metadata = pd.Series(name=json.dumps(metadata), dtype='object')
seqsdata = pd.DataFrame([["EVQLVESGGGLAKPGGSLRLHCAASGFAFSSYWMNWVRQAPGKRLEWVSAINLGGGLTYYAASVKGRFTISRDNSKNTLSLQMNSLRAEDTAVYYCATDYCSSTYCSPVGDYWGQGVLVTVSS"],
                          ["EVQLVQSGAEVKRPGESLKISCKTSGYSFTSYWISWVRQMPGKGLEWMGAIDPSDSDTRYNPSFQGQVTISADKSISTAYLQWSRLKASDTATYYCAIKKYCTGSGCRRWYFDLWGPGT"]
                         ], columns = ['sequence_alignment_aa'])
                         
save_file = "../data/custom-data-example-2.csv"
metadata.to_csv(save_file, index=False)
seqsdata.to_csv(save_file, index=False, mode='a')

```

### 2. Turn your OAS formatted files into a custom database

After creating all the files you want to include in the new database, you can run the following code to create the database.

```{r, engine='python', count_lines}
db_folder = "../data/my_db"
db_files = ['../data/custom-data-example.csv']

newDB = PrepareDB(db_path=db_folder)

db_dict = {}
for num, data_unit_file in enumerate(db_files):

    db_dict[num] = data_unit_file
    
    metadata = json.loads(','.join(pd.read_csv(data_unit_file, nrows=0).columns))
    seqsdata = pd.read_csv(data_unit_file, header=1, usecols=['sequence_alignment_aa']).iloc[:,0].values
    numbered_sequences = [number(sequence, allowed_species=None)[0] for sequence in seqsdata]

    newDB.prepare_database(numbered_sequences, 
                           file_id=num, 
                           chain=metadata['Chain'], 
                           species=metadata['Species'])
    
newDB.save_database()

merge_files(db_folder)

with open(os.path.join(db_folder, "my_db_id_to_study.txt"), "w") as handle: handle.write(str(db_dict))

```

### 3. Initiate the search class with your custom database

```{r, engine='python', count_lines}
mydb = SearchDB(database_path=db_folder,      # Path to your database. Default will be to download a prepared version of OAS.
                 allowed_chain='Heavy',       # Search against a specific chain. Default is any chain.
                 allowed_species='Any',       # Search against a specific species. Default is any species.
                 n_jobs=1                     # Allocated number for jobs/threads for the search.
                )
```


-----



### Citation   
Work in preparation
