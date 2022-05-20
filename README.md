
---

<div align="center">    
 
# KA-Search: Rapid and exact sequence identity search of known antibodies  


</div>


In antibody therapeutic drug discovery, finding antibodies with similar sequence, or complementary determining regions (CDR) as an antibody of interest is a convenient method for understanding ones drug. A simple and often used approach for finding similar sequences is using sequence identity, either for the whole sequence or CDRs. However, with the current and ever-growing size of available antibody sequences, searching against all known antibodies is becoming ever the more difficult. The Observed Antibody Space (OAS) database currently contains 1.7 billion antibody sequences. A tool optimized for this specific task is therefore relevant for speeding up this process and making it feasible.

Here, we introduce Known Antibody Search (KA-Search), a tool for rapid search of similar whole chain, CDRs and CDR3 antibodies in the OAS database. We demonstrate how with a prepared target dataset, KA-Search can be used to find the N most similar antibodies from a dataset of 1.7 billion antibodies within 1.5 hour on a \_ CPU. Further, for convenience, a user can create their own database to search against, utilizing the same functionality and speed on in-house data. 

-----------

# Install KA-Search

KA-Search is freely available and can be installed with pip.

~~~.sh
    pip install kasearch
~~~

or directly from github.

~~~.sh
    pip install -U git+
~~~

----------

# Searching with KA-Search

**A Jupyter notebook** showcasing KA-Search can be found [here](https://github.com/TobiasHeOl/AbLang/tree/main/examples). 

First, your sequence(s) needs to be transformed into the correct alignment.

```{r, engine='python', count_lines}

from kasearch import SetQueries

raw_queries = [
    'QVKLQESGAELARPGASVKLSCKASGYTFTNYWMQWVKQRPGQGLDWIGAIYPGDGNTRYTHKFKGKATLTADKSSSTAYMQLSSLASEDSGVYYCARGEGNYAWFAYWGQGTTVTVSS',
    'EVQLQQSGTVLARPGASVKMSCEASGYTFTNYWMHWVKQRPGQGLEWIGAIYPGNSDTSYIQKFKGKAKLTAVTSTTSVYMELSSLTNEDSAVYYCTLYDGYYVFAYWGQGTLVTVSA',
    'QVQLLESGGGLVQPGGSLRLSCAASGFTFSTAAMSWVRQAPGKGLEWVSGISGSGSSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARELSYLYSGYYFDYWGQGTLVTVSS',
]

QueryDB = SetQueries(raw_queries)
QueryDB.queries[0]
```
Acheiving an output that looks like this:

```console
Query(aligned_query=array([81, 86, 75,  0, 76, 81, 69, 83, 71, 65,  0, 69, 76, 65, 82, 80, 71,
       65, 83, 86, 75, 76, 83, 67, 75, 65, 83, 71, 89, 84, 70,  0,  0,  0,
        0,  0,  0,  0,  0,  0, 84, 78, 89, 87, 77, 81,  0, 87, 86, 75, 81,
        0, 82,  0, 80,  0, 71,  0, 81,  0,  0, 71,  0, 76, 68,  0, 87, 73,
       71, 65, 73, 89, 80, 71,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       68, 71, 78, 84, 82, 89,  0,  0, 84,  0,  0, 72,  0,  0, 75, 70,  0,
        0, 75,  0,  0,  0, 71, 75, 65, 84, 76, 84, 65,  0, 68,  0,  0,  0,
       75,  0, 83,  0,  0, 83, 83,  0,  0,  0,  0, 84,  0, 65, 89, 77, 81,
       76, 83, 83, 76, 65, 83,  0, 69, 68, 83, 71, 86, 89, 89, 67, 65, 82,
       71, 69, 71, 78,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 89, 65, 87, 70, 65,
       89, 87, 71,  0, 81, 71, 84, 84, 86, 84, 86, 83, 83], dtype=int8), chain='Heavy')
```

You then initiate your search and search against the target database of choice.

```{r, engine='python', count_lines}
search_oas = SearchOAS(database_path=DATABASE_PATH, 
                       allowed_chain='Heavy', 
                       n_jobs=2
                      )
                      
search_oas.search(QueryDB.queries[0], keep_best_n=2)
search_oas.current_best_identities
```

Returning the following best identities

```console
array([[0.79510413, 0.78571429, 0.83333333],
       [0.79161837, 0.78571429, 0.78409091]])
```

You can then access the indices and use them to extract the meta data of the N best matches.


-----



### Citation   
Work in preparation
