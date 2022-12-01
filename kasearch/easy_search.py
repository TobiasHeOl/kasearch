from kasearch import AlignSequences, SearchDB


def EasySearch(query, 
               keep_best_n=10,
               database_path='oasdb-tiny', 
               allowed_chain='Any', 
               allowed_species='Any',
               regions=['whole'],
               length_matched=[False],
               n_jobs=1,
              ):
    """Quick KA-Search wrapper to run of a single query across a single region. 

    Parameters
    ----------
    query : str
        Antibody sequence to search with
    keep_best_n : int
        Number of closest matches to return (default is 10)
    database_path : str
        Path to processed database to search against (default is oasdb-small)
    allowed_chain : str
        Which chain to search against (default is Any, options are Any, Heavy or Light)
    allowed_species : str
        Which species to search against (default is Any, options depend on database searched but often include Human and Mouse)
    regions : list
        List of regions to search. Only the first will be used (default is [whole])
    length_matched : list
        A flag for whether the identity will only be calculated between regions of identical length (default is [False])
    n_jobs : int
        Number of threads used (default is 1)

    Returns
    -------
    pandas dataframe
        A pandas dataframe with the keep_best_n closest sequences, and their meta data, to the query
    """
    
    querydb = AlignSequences(n_jobs=n_jobs)(query)
    
    targetdb = SearchDB(database_path=database_path,
                    allowed_chain=allowed_chain, 
                    allowed_species=allowed_species, 
                    regions=regions, 
                    length_matched=length_matched,
                   )
    
    targetdb.search(querydb[:1], keep_best_n=keep_best_n)

    return targetdb.get_meta(n_query=0, n_region=0, n_sequences='all', n_jobs=n_jobs)