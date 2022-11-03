import re
import numpy as np

# This has to be a numbering large enough to fit any sequence you will use. If a given number is not here it will break
canonical_numbering = ["1 ", "2 ", 
                       "3 ", "3A", 
                       "4 ", "5 ", "6 ", "7 ", "8 ", "9 ", "10 ", "11 ", "12 ", "13 ", "14 ", "15 ", "16 ", "17 ", "18 ", 
                       "19 ", "20 ", "21 ", "22 ", "23 ", "24 ", "25 ", "26 ", "27 ", "28 ", "29 ", "30 ", "31 ",  
                       "32 ", "32A", "32B", '33C', "33B", "33A", "33 ", 
                       "34 ", "35 ", "36 ", "37 ", "38 ", "39 ", 
                       "40 ", '40A', 
                       "41 ", "42 ", "43 ", 
                       "44 ", "44A",
                       "45 ", "45A", 
                       "46 ", "46A",
                       "47 ", "47A", 
                       "48 ", "48A", "48B", 
                       "49 ", "49A", 
                       "50 ", 
                       "51 ", "51A", 
                       "52 ", "53 ", "54 ", "55 ", "56 ", "57 ", "58 ", "59 ", 
                       "60 ", "60A", "60B", "60C", '60D', "61E", '61D', "61C", "61B", "61A", "61 ", 
                       "62 ", "63 ", "64 ", "65 ", "66 ", 
                       "67 ", '67A', "67B", 
                       "68 ", "68A", "68B", 
                       "69 ", "69A", "69B",
                       "70 ", 
                       "71 ", "71A", "71B", 
                       "72 ", 
                       "73 ", '73A', "73B",
                       "74 ", "75 ", "76 ", "77 ", "78 ", "79 ", 
                       "80 ", "80A", 
                       "81 ", "81A", "81B", "81C",
                       "82 ", "82A", 
                       "83 ", "83A", "83B",
                       "84 ", 
                       "85 ", "85A", "85B", "85C", "85D", 
                       "86 ", "86A", 
                       "87 ", "88 ", "89 ", "90 ", "91 ", "92 ", "93 ", "94 ", "95 ", 
                       "96 ", "96A",
                       "97 ", "98 ", "99 ", "100 ", "101 ", "102 ", "103 ", "104 ", "105 ", "106 ", "107 ", "108 ", 
                       "109 ", "110 ",
                       "111 ", "111A", "111B", "111C", "111D", "111E", "111F", "111G", '111H', '111I', '111J', 
                       "111K", "111L", "112L",
                       '112K', '112J', '112I', '112H', "112G", "112F", "112E", "112D", "112C", "112B", "112A", "112 ",
                       "113 ","114 ","115 ","116 ","117 ","118 ",
                       "119 ", "119A",
                       "120 ","121 ","122 ","123 ","124 ","125 ", "126 ","127 ","128 "
                        ]

canonical_numbering_to_index = {x: i for i, x in enumerate(canonical_numbering)}
canonical_numbering_len = len(canonical_numbering)

# SEQUENCES ARE NUMBERED USING IMGT NUMBERING
reg_def = dict()
reg_def["CDR1"] = list(range(27, 38 + 1))
reg_def["CDR2"] = list(range(56, 65 + 1))
reg_def["CDR3"] = list(range(105, 117 + 1))
reg_def["CDR_all"] = reg_def["CDR1"] + reg_def["CDR2"] + reg_def["CDR3"]

cdr3_mask = np.array([int(canonical_numbering[i][:-1]) in reg_def["CDR3"] for i in range(len(canonical_numbering))], dtype = np.uint8)
all_cdrs_mask = np.array(
    [int(canonical_numbering[i][:-1]) in reg_def["CDR_all"] for i in range(len(canonical_numbering))], dtype = np.uint8)


def get_region_mask(region):
    """Retrieve mask for a specific region.
    
    Parameters
    ----------
    region : list
        List of unique positions from the canonical alignment to not mask

    Returns
    -------
    numpy array
        Specific region mask
    """     
    
    if region == 'whole':
        return np.ones(200, dtype = np.uint8)
    elif region == 'cdrs':
        return all_cdrs_mask
    elif region == 'cdr3':
        return cdr3_mask
    else:
        assert np.all([i in canonical_numbering for i in region]), "At least one position in the user defined region, is not part of the 200 used positions. \
        Please remove such positions."
        return np.array([i in region for i in canonical_numbering], dtype = np.uint8)


def canonical_alignment(anarci_output):
    """Alignment of ANARCI output.
    
    Parameters
    ----------
    anarci_output : dict
        ANARCI numbering of a sequence

    Returns
    -------
    numpy array
        Canonical alignment of a sequence
    """

    sequence = np.zeros(canonical_numbering_len, np.int8)

    for res in anarci_output:
        num, seq = res
        num = ''.join([str(i) for i in num])

        i = canonical_numbering_to_index[num]
        sequence[i] = int.from_bytes(seq.encode(), 'big')

    sequence[(sequence == 45)] = 0
    return sequence


oas_numbering_finder = re.compile('\d+.\'\: \'[A-Z]')
def canonical_alignment_oas(anarci_output):
    """OAS derived ANARCI output.
    
    Parameters
    ----------
    anarci_output : dict
        OAS derived ANARCI numbering of a sequence

    Returns
    -------
    numpy array
        Canonical alignment of a sequence
    """
    
    sequence = np.zeros(canonical_numbering_len, np.int8)

    for res in oas_numbering_finder.findall(anarci_output):
        num, seq = res[:-5], res[-1]
        i = canonical_numbering_to_index[num]
        sequence[i] = int.from_bytes(seq.encode(), 'big')

    sequence[(sequence == 45)] = 0
    return sequence