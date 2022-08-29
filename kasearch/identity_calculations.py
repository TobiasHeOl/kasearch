import jax
import jax.numpy as jnp

import numpy as np 
from multiprocessing import Pool
from kasearch.canonical_alignment import all_cdrs_mask, cdr3_mask, reg_def

default_region_masks = np.stack([np.ones(200, dtype = np.uint8), all_cdrs_mask, cdr3_mask])
default_length_matched = np.array([False,True,True], dtype = bool)

@jax.jit
def calculate_seq_ids_multiquery(array_of_abs1, array_of_abs2, region_masks, length_matched):

    masks, length_matched = jnp.array(region_masks.T), jnp.array(length_matched)
    
    abs1 = jax.lax.stop_gradient(jnp.expand_dims(array_of_abs1, axis=1))
    abs2 = jax.lax.stop_gradient(jnp.expand_dims(array_of_abs2, axis=0))
    comparison = abs1 == abs2

    mask1, mask2 = abs1 != 0, abs2 != 0
    overlapping_residues = comparison * mask1 * mask2

    len1, len2 = mask1 @ masks, mask2 @ masks
    region_overlap = overlapping_residues @ masks
    identities = ((region_overlap / len1) + (region_overlap / len2)) * (~length_matched | (len1==len2)) / 2

    return identities



def get_n_most_identical_multiquery(query, targets, target_ids, n=10,
                                    region_masks=default_region_masks, 
                                    length_matched=default_length_matched):
    
    n = len(targets)-1 if len(targets) < n else n # Adjusts for large n's
    
    seq_identity_matrix = calculate_seq_ids_multiquery(query, targets, region_masks, length_matched)
    seq_identity_matrix = np.array(seq_identity_matrix)
    seq_identity_matrix[np.isnan(seq_identity_matrix)] = 0

    position_of_n_best = np.argpartition(-seq_identity_matrix, n, axis=1)  # partition by seq_id
    n_highest_identities = np.take_along_axis(seq_identity_matrix, position_of_n_best, axis=1)[:,:n]

    broadcasted_ids = np.broadcast_to(target_ids[:, None], (query.shape[0], targets.shape[0], region_masks.shape[0], 2))

    n_highest_ids = np.take_along_axis(broadcasted_ids, position_of_n_best[:, :, :, None], axis=1)[:,:n]

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
    n_jobs = n_jobs if n_jobs is not None else numba.get_num_threads() - 1
    seq_identity_matrix = np.array(slow_calculate_many_seq_ids(query, target, n_jobs=n_jobs))

    position_of_n_best = np.argpartition(-seq_identity_matrix, n, axis=0)  # partition by seq_id
    n_highest_identities = np.take_along_axis(seq_identity_matrix, position_of_n_best, axis=0)[:n]

    broadcasted_ids = np.broadcast_to(target_ids[:, None], (len(target), 3, 2))
    n_highest_ids = np.take_along_axis(broadcasted_ids, position_of_n_best[:, :, None], axis=0)[:n]

    return n_highest_identities, n_highest_ids
