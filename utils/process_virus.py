from utils.virus_structure import *
from utils.utils import *
from collections import Counter
import pandas as pd


def get_virus(filename):
    """
    Reads the given virus filename and processes it
    
    Parameters:
        filename (String): fasta file
    """
    virus = read_file(filename)
    (virus.codon_freq, virus.decodon_freq) = get_virus_info(virus)
    return virus


def get_frames_start_stop_pairs(virus_sequence, name, reverse_seq = False, print_results = False):
    """
    Splits the given virus sequences into frames and calculates all start stop pairs without stop codon in between
    
    Parameters:
        virus_sequence (Seq): original virus sequence to split into frames
        name (Dict): set of frames
        reverse_seq (Boolean): print results for the reversed frame 
        print_results (Boolean): print results into console
    """
    i = 0
    for key in name:
        name[key]['seq'] = virus_sequence[i:]
        name[key]['codons'] = get_codons(name[key]['seq'])
        start_stop_pairs = get_start_stop_pairs(name[key]['codons'])
        name[key]['start_stop_p'] = start_stop_pairs
        
        if print_results:
            if reverse_seq:
                print(f'Seka: {reverse}; skaitymo remelis: {key[-1]}; Start/stop kodonų poros: {start_stop_pairs}\n') 

            else:
                print(f'Seka: {straight}; skaitymo remelis: {key[-1]}; Start/stop kodonų poros: {start_stop_pairs}\n') 

        i += 1


def filter_long_fragments(virus_frames, reverse_seq = False):
    """
    Returns fragments longer than 100 symbols with cropped stop codon for the set of frames
    
    Parameters:
        virus_frames (Dict): dictionary of frames
    """
    for key in virus_frames:
        long_frag = get_long_fragments(virus_frames[key]['seq'], virus_frames[key]['start_stop_p'])
        if reverse_seq:
            print(f'Seka: {reverse}; skaitymo remelis: {key[-1]}; Fragmentai ilgesni nei 100 simboliu: {long_frag}\n') 

        else:
            print(f'Seka: {straight}; skaitymo remelis: {key[-1]}; Fragmentai ilgesni nei 100 simboliu: {long_frag}\n') 


def get_frames_codon_freq(virus_frames):
    """
    Calculates codon freq for the set of frames
    
    Parameters:
        virus_frames (Dict): dictionary of frames
    """
    for key in virus_frames:
        codon_freq = Counter(virus_frames[key]['codons'])
        l = seq_codons[key]
        for el in l:
            l[el] = codon_freq[el]


def get_frames_decodon_freq(virus_frames):
    """
    Calculates decodon freq for the set of frames
    
    Parameters:
        virus_frames (Dict): dictionary of frames
    """
    for key in virus_frames:
        virus_frames[key]['decodon_freq'] = Counter(get_decodons(virus_frames[key]['seq']))
        l = seq_decodons[key]
        for el in l:
            l[el] = virus_frames[key]['decodon_freq'][el]


def get_total_codon_freq(seq_codons, show_table = False):
    """
    Returns the total virus codon frequency table
    
    Parameters:
        seq_codons (Dict): dictionary of frames with their respective codon frequency
        show_table (Boolean): print results into console
    """
    total_codon_freq = Counter()
    for frame in seq_codons:
        total_codon_freq.update(seq_codons[frame])

    for (frame, codon_dict) in seq_codons.items():
        for c in codon_dict:
            if total_codon_freq[c] == 0:
                codon_dict[c] = 0
            else:
                codon_dict[c] = codon_dict[c] / total_codon_freq[c]

    if show_table:
        df = pd.DataFrame(seq_codons)
        print(f'Kodonu daznis: \n\n{df}')

    return total_codon_freq


def get_total_decodon_freq(seq_decodons, show_table = False):
    """
    Returns the total virus decodon frequency table
    
    Parameters:
        seq_decodons (Dict): dictionary of frames with their respective decodon frequency
        show_table (Boolean): print results into console
    """
    total_decodon_freq = Counter()
    for frame in seq_decodons:
        total_decodon_freq.update(seq_decodons[frame])
        
    for (frame, decodon_dict) in seq_decodons.items():
        for c in decodon_dict:
            if total_decodon_freq[c] == 0:
                decodon_dict[c] = 0
            else:
                decodon_dict[c] = decodon_dict[c] / total_decodon_freq[c]

    if show_table:
        df = pd.DataFrame(seq_decodons)
        print(f'Dikodonu daznis: \n\n{df}')

    return total_decodon_freq


def get_virus_info(self):
    """
    Returns the information of the virus needed for multiple virus comparison
    
    Parameters:
        self (Seq): original virus sequence
    """
    virus_seq = self.seq
    reverse_virus_sequence = reverse_sequence(virus_seq)
    virus_frames = seq_frames
    reverse_virus_frames = reverse_seq_frames

    get_frames_start_stop_pairs(virus_seq, virus_frames)
    get_frames_start_stop_pairs(reverse_virus_sequence, reverse_virus_frames)

    get_frames_codon_freq(virus_frames)
    get_frames_codon_freq(reverse_virus_frames)

    get_frames_decodon_freq(virus_frames)
    get_frames_decodon_freq(virus_frames)

    total_codon_freq = get_total_codon_freq(seq_codons)
    total_decodon_freq = get_total_decodon_freq(seq_decodons)
    return (total_codon_freq, total_decodon_freq)