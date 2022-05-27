import numpy as np
import numba
from multiprocessing import Pool
from kasearch.canonical_alignment import all_cdrs_mask, cdr3_mask, reg_def


@numba.njit("UniTuple(float32, 3)(int8[:],int8[:])",error_model="numpy", fastmath=True, cache=True)
def calculate_seq_id(ab1, ab2):
    """ Takes two canonically aligned sequences and computes sequence identity
    Only used for non-parallel.
    """
    comparison = ab1 == ab2
    mask1, mask2 = ab1 != 0, ab2 != 0
    overlapping_residues = comparison * mask1 * mask2

    # For the whole sequence:
    full_overlap = np.sum(overlapping_residues)
    full_id = (full_overlap / mask1.sum() + full_overlap / mask2.sum()) / 2
        
    # For all CDRs
    cdrs_overlap = np.sum(all_cdrs_mask * overlapping_residues)
    cdrs_id = (cdrs_overlap / (mask1 * all_cdrs_mask).sum() + cdrs_overlap / (mask2 * all_cdrs_mask).sum()) / 2
    
    # For CDR-H3
    h3_len1, h3_len2 = (mask1 * cdr3_mask).sum(), (mask2 * cdr3_mask).sum()
    h3_overlap = np.sum(cdr3_mask * overlapping_residues)
    h3_id = (h3_overlap / h3_len1 + h3_overlap / h3_len2) / 2
   
    return full_id, cdrs_id, h3_id


@numba.njit("float32[:,:](int8[:],int8[:,:])",error_model="numpy", fastmath=True, cache=True)
def calculate_all_seq_ids(ab1, array_of_abs):
    """ Computes sequence identity of one antibody sequence against an array of sequences
    Only used for non-parallel.
    """
    
    size = array_of_abs.shape[0]

    identities = np.empty((size, 3), dtype=np.float32)   
    
    for i in range(size):
        identities[i] = calculate_seq_id(ab1, array_of_abs[i])

    return identities


@numba.njit("float32[:,:](int8[:],int8[:,:],float32[:,:])",error_model="numpy", fastmath=True, cache=True)
def _calculate_batch_seq_ids(ab1, ab2, identities):
    """ Computes sequence identity of one antibody sequence against an array of sequences
    Not meant to be used outside of calculate_all_seq_ids_parallel
    """
    comparison = ab1 == ab2
    mask1, mask2 = ab1 != 0, ab2 != 0
    overlapping_residues = comparison * mask1 * mask2

    # For the whole sequence:
    full_overlap = np.sum(overlapping_residues, axis=-1)
    identities[:,0] = (full_overlap / mask1.sum() + full_overlap / mask2.sum(axis=-1)) / 2
        
    # For all CDRs
    cdrs_overlap = np.sum(all_cdrs_mask * overlapping_residues, axis=-1)
    identities[:,1] = (cdrs_overlap / (mask1 * all_cdrs_mask).sum() + cdrs_overlap / (mask2 * all_cdrs_mask).sum(axis=-1)) / 2
    
    # For CDR-H3
    h3_len1, h3_len2 = (mask1 * cdr3_mask).sum(), (mask2 * cdr3_mask).sum(axis=-1)
    h3_overlap = np.sum(cdr3_mask * overlapping_residues, axis=-1)
    identities[:,2] = ((h3_overlap / h3_len1 + h3_overlap / h3_len2) / 2) * (h3_len1==h3_len2)
   
    return identities


@numba.njit("float32[:,:](int8[:],int8[:,:])",error_model="numpy", fastmath=True, parallel=True)
def calculate_all_seq_ids_parallel(ab1, array_of_abs):
    """ Computes sequence identity of one antibody sequence against an array of sequences
    Only used for parallel.
    """
    size = array_of_abs.shape[0]
    chunk_size = 200  # Somewhat arbitrary number (too big is slow and too small is slow)
    chunks = size//chunk_size
    
    identities = np.empty((size, 3), dtype = np.float32)
    global_buffer = np.empty((numba.get_num_threads(), chunk_size, 3), dtype = np.float32)
    
    for i in numba.prange(chunks):
        thread_id = numba.np.ufunc.parallel._get_thread_id()
        start, end = i*chunk_size, (i+1)*chunk_size
        thread_buffer = global_buffer[thread_id, :size-end]
        identities[start:end] = _calculate_batch_seq_ids(ab1, array_of_abs[start:end], thread_buffer)

    return identities


def get_n_most_identical(query, target, target_ids, n=10, n_jobs=None):
    if n_jobs != 1:
        if n_jobs is not None:
            numba.set_num_threads(n_jobs)
        seq_identity_matrix = calculate_all_seq_ids_parallel(query, target)
    else:
        seq_identity_matrix = calculate_all_seq_ids(query, target)
        
    where_are_NaNs = np.isnan(seq_identity_matrix)
    seq_identity_matrix[where_are_NaNs] = 0

    position_of_n_best = np.argpartition(-seq_identity_matrix, n, axis=0)  # partition by seq_id
    n_highest_identities = np.take_along_axis(seq_identity_matrix, position_of_n_best, axis=0)[:n]

    broadcasted_ids = np.broadcast_to(target_ids[:, None], (target.shape[0], 3, 2))
    n_highest_ids = np.take_along_axis(broadcasted_ids, position_of_n_best[:, :, None], axis=0)[:n]

    return n_highest_identities, n_highest_ids


def slow_calculate_seq_id(ab1, ab2):
    # Around 1k times slower but doesn't rely on canonical alignment
    ab1_clean = [x for x in ab1 if x[1] != "-"]
    ab2_clean = [x for x in ab2 if x[1] != "-"]

    overlapping_residues = [x for x in ab1_clean if x in ab2_clean]

    # For the whole sequence:
    full_overlap = len(overlapping_residues)
    full_id = (full_overlap / len(ab1_clean) + full_overlap / len(ab2_clean)) / 2

    # For all CDRs
    cdrs_overlap = len([x for x in overlapping_residues if x[0][0] in reg_def["CDR_all"]])
    cdrs_len1 = len([x for x in ab1_clean if x[0][0] in reg_def["CDR_all"]])
    cdrs_len2 = len([x for x in ab2_clean if x[0][0] in reg_def["CDR_all"]])

    if cdrs_len1 and cdrs_len2:
        cdrs_id = (cdrs_overlap / cdrs_len1 + cdrs_overlap / cdrs_len2) / 2
    else:
        cdrs_id = 0

    # For CDR-H3
    h3_overlap = len([x for x in overlapping_residues if x[0][0] in reg_def["CDR3"]])
    h3_len1 = len([x for x in ab1_clean if x[0][0] in reg_def["CDR3"]])
    h3_len2 = len([x for x in ab2_clean if x[0][0] in reg_def["CDR3"]])

    if h3_len1 and h3_len2:
        h3_id = (h3_overlap / h3_len1 + h3_overlap / h3_len2) / 2
    else:
        h3_id = 0

    return full_id, cdrs_id, h3_id


def slow_calculate_many_seq_ids(ab1, list_of_abs, n_jobs=1):
    size = len(list_of_abs)

    with Pool(processes=n_jobs) as pool:
        return pool.starmap(slow_calculate_seq_id, zip(size * [ab1], list_of_abs), chunksize=size // n_jobs)


def slow_get_n_most_identical(query, target, target_ids, n=10, n_jobs=None):
    n_jobs = n_jobs if n_jobs is not None else numba.get_num_threads()
    seq_identity_matrix = np.array(slow_calculate_many_seq_ids(query, target, n_jobs=n_jobs))

    position_of_n_best = np.argpartition(-seq_identity_matrix, n, axis=0)  # partition by seq_id
    n_highest_identities = np.take_along_axis(seq_identity_matrix, position_of_n_best, axis=0)[:n]

    broadcasted_ids = np.broadcast_to(target_ids[:, None], (len(target), 3, 2))
    n_highest_ids = np.take_along_axis(broadcasted_ids, position_of_n_best[:, :, None], axis=0)[:n]

    return n_highest_identities, n_highest_ids
