from kasearch import AlignSequences, SearchDB


def EasySearch(query, 
               keep_best_n=10,
               database_path='oasdb-small', 
               allowed_chain='Any', 
               allowed_species='Any',
               regions=['whole'],
               length_matched=[False],
               n_jobs=1,
              ):
    
    querydb = AlignSequences(n_jobs=n_jobs)(query)
    
    targetdb = SearchDB(database_path=database_path,
                    allowed_chain=allowed_chain, 
                    allowed_species=allowed_species, 
                    regions=regions, 
                    length_matched=length_matched,
                   )
    
    targetdb.search(querydb[:1], keep_best_n=keep_best_n)

    return targetdb.get_meta(n_query=0, n_region=0, n_sequences='all', n_jobs=n_jobs)