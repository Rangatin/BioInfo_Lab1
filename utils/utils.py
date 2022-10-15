import os
import pandas as pd
from Bio import SeqIO
from itertools import product
from scipy.spatial.distance import pdist, squareform

symbols = 'ACGT'
start_codon = 'ATG'
stop_codons = ['TGA', 'TAG', 'TAA']

straight = 'tiesiogine'
reverse = 'reverse komplementas'

# file reading config
__location__ = os.path.realpath(os.getcwd())


def read_file(filename):
    """
    Returns the parsed sequence
    
    Parameters:
        filename (String): name of the file to parse
    """
    with open(os.path.join(__location__, 'data', filename)) as handle:
        return SeqIO.read(handle, 'fasta')
    

def reverse_sequence(seq):
    """
    Returns the reverse complement of a sequence
    
    Parameters:
        seq (Seq): a sequence
    """
    return seq.reverse_complement()


def get_all_possible_codons():
    """
    Returns: frequency dictionary of all possible codons
    """
    return {''.join(k): 0 for k in product(symbols, repeat=3)} 


def get_all_possible_decodons():
    """
    Returns: frequency dictionary of all possible decodons
    """
    return {''.join(k): 0 for k in product(symbols, repeat=6)} 


def get_codons(seq):
    """
    Returns: given sequence split into codons
    
    Parameters:
        seq (Seq): a sequence
    """
    n = 3 # codons
    found_codons = []
    for i in range(0, len(seq) - n + 1, n):
        found_codons.append(str(seq[i:i+n]))
    
    return found_codons


def get_decodons(seq, reverse = False):
    """
    Returns: given sequence split into decodons
    
    Parameters:
        seq (Seq): a sequence
    """
    n = 6 # decodons
    found_codons = []
    if reverse:
        for i in reversed(range(0, len(seq) - n + 1, n)):
                found_codons.append(str(seq[i:i+n]))
    
    else:
        for i in range(0, len(seq) - n + 1, n):
            found_codons.append(str(seq[i:i+n]))
    
    return found_codons

def get_start_stop_pairs(codons):
    """
    Returns: tuple indexes of start and stop codon pairs
    
    Parameters:
        codons (List): a list of codons
    """
    if(len(codons) == 0):
        return []

    else:
        start_codon_present = False
        found_pairs = []

        for idx in range(0, len(codons)):
            if(not start_codon_present):
                if(codons[idx] == start_codon):
                    start_idx = idx * 3
                    start_codon_present = True

            elif(codons[idx] in stop_codons):
                found_pairs.append((start_idx, idx * 3 + 2))
                start_codon_present = False

        return found_pairs


def get_long_fragments(seq, start_stop_pairs):
    """
    Returns: sequence fragments longer than 100 symbols
    
    Parameters:
        seq (Seq): a sequence
        start_stop_pairs(List): tuple indexes of start and stop codon pairs
    """
    if len(start_stop_pairs) == 0:
        return []

    else:
        long_fragments = []
        for (start_idx, stop_idx) in start_stop_pairs:
            if(stop_idx - start_idx + 1 > 100):
                long_fragments.append(seq[start_idx:stop_idx])
                
        return long_fragments

        
def get_distance_matrix(freq_table):
    """
    Converts a given frequency table into a distance matrix
    
    Parameters:
        freq_table (Dict): frequency table
    """
    virus_df = pd.DataFrame(freq_table).T
    pairwise = pd.DataFrame(
        squareform(pdist(virus_df)),
        columns = virus_df.index,
        index = virus_df.index
    )
    return pairwise


def write_to_phylip(filename, dist_matrix):
    """
    Writes a given distance matrix to a phylip file
    
    Parameters:
        filename (String): phylip file
        dist_matrix(DataFrame): distance matrix to write
    """
    with open(os.path.join('res', filename), 'w') as f:
        f.write('%d\n' % len(dist_matrix))
        dfAsString = dist_matrix.to_string(header=False)
        f.write(dfAsString)

