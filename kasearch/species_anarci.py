from anarci import anarci
from multiprocessing import Pool
from functools import partial

def _number_chunk(sequences, scheme="imgt", database="ALL", allow=set(["H","K","L"]), allowed_species=['human','mouse'], strict = True, **kwargs):
    try:
        numbered, _, _ = anarci(list(enumerate(sequences)), scheme=scheme, database=database, allow=allow, allowed_species=allowed_species, **kwargs)
        numbered  = [x[0][0] if x else None for x in numbered]
    except Exception:
        numbered = []
        for seq in sequences:
            try:
                numb, _, _ = anarci([("sequence",seq)], scheme=scheme, database=database, allow=allow, allowed_species=allowed_species,**kwargs)

                if numb[0]:
                    numb = numb[0][0][0]
                else:
                    numb = None
            except Exception:
                numb = None

            numbered.append(numb)

    if strict:
        assert None not in numbered, "Failed numbering for at least one sequence in the set"

    return numbered

def number_many_at_once(sequences, scheme="imgt", database="ALL", allow=set(["H","K","L"]), allowed_species=['human','mouse'], strict=True, ncpu=1, **kwargs):
    chunk_size = 10_000
    kwargs['ncpu'] = 1 #For HMMER 
    kwargs['output'] = False

    number_chunk = partial(_number_chunk, scheme=scheme, database=database, allow=allow, allowed_species=allowed_species, strict=strict, **kwargs)

    if ncpu==1:
        numbered = sum([number_chunk(sequences[i:i+chunk_size]) for i in range(0, len(sequences), chunk_size)], [])
    else:
        with Pool(ncpu) as pool:
            numbered = pool.map(number_chunk, [sequences[i:i+chunk_size] for i in range(0, len(sequences), chunk_size)], chunksize=1)
    
    return numbered
            
    
