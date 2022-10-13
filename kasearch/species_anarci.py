from anarci import validate_sequence, anarci, chain_type_to_class, scheme_short_to_long, run_anarci


# Wrapper function for simple sequence in numbering and chain type out behaviour.
def number(sequence, scheme="imgt", database="ALL", allow=set(["H", "K", "L"]), allowed_species=None):
    """
    Given a sequence string, use anarci to number it using the scheme of choice.
    Only the first domain will be recognised and numbered
    For multiple sequences it is advised to use run_anarci instead of iterative use of this function.
    @param sequence: An amino acid sequence string
    @param scheme: The numbering scheme that should be applied. Choose from imgt, chothia, kabat or martin
    @param database: The HMMER database that should be used. Normally not changed unless a custom db is created.
    @param allow: A set containing the chain types that should be recognised. If chothia, kabat or martin is used
                  as the scheme, anarci will ignore tcr chains.
    @param allowed_species: A set containing the species that should be recognised.
    @return: If the sequence can be numbered, a list containing the numbering and sequence; and the chain type.
             Otherwise both are False.
    """

    try:
        validate_sequence(sequence)
        scheme = scheme_short_to_long[scheme.lower()]
    except KeyError:
        raise AssertionError("Unrecognised to unimplemented scheme: %s" % scheme)

    if len(
            sequence) < 70:
        return False, False

    try:
        numbered, alignment_details, _ = anarci([("sequence_0", sequence)], scheme=scheme, database=database,
                                                output=False, allow=allow, allowed_species=allowed_species)
    except AssertionError:  # Catch where the user has tried to number a TCR with an antibody scheme
        return False, False

    # We return the numbering list and the chain type where kappa and lambda chains are both "L" for light
    if numbered[0]:
        return numbered[0][0][0], get_chain(alignment_details)
    else:
        return False, False

def get_chain(alignment_details):
    
    chain = chain_type_to_class[alignment_details[0][0]["chain_type"]]
    
    if chain == 'H':
        return "Heavy"
    
    elif (chain == 'L') or (chain == 'K'):
        return "Light"
    
    
def many_number(sequences, scheme="imgt", database="ALL", allow=set(["H", "K", "L"]), allowed_species=None, n_jobs=1):
    """
    Less robust, but much faster anarci numbering.
    """
    sequences = [("sequence_{}".format(num), sequence) for num, sequence in enumerate(sequences)]   
    
    try:
        _, numbered_seqs, _, _ = run_anarci(sequences, 
                                     scheme=scheme, 
                                     database=database,
                                     allow=allow, 
                                     allowed_species=allowed_species,
                                     ncpu=n_jobs,
                                    )
        
        return numbered_seqs
    
    except Exception:  # Catch where the user has tried to number a TCR with an antibody scheme
        return False
            
    