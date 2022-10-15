from utils.utils import *
from utils.virus_structure import *
from utils.process_virus import *
import sys

def main():
    """
    1. Pateiktoje sekoje fasta formatu surastu visas start ir stop kodonų poras,
    tarp kurių nebutu stop kodono (ir tiesioginei sekai ir jos reverse komplementui). 
    2. Kiekvienam stop kodonui parinkti toliausiai nuo jo esanti start kodoną (su salyga, kad tarp ju nera kito stop kodono)
    """

    filename = sys.argv[1]
    virus = read_file(filename)

    virus_sequence = virus.seq
    reverse_virus_sequence = reverse_sequence(virus_sequence)

    get_frames_start_stop_pairs(virus_sequence, seq_frames, print_results = True)
    print('\n\n')
    get_frames_start_stop_pairs(reverse_virus_sequence, reverse_seq_frames, reverse_seq = True, print_results = True)
    print('\n\n')

    """
    3. Atfiltruokite visus fragmentus ("tai butu baltymų koduojancios sekos"), kurie trumpesni nei 100 fragmentų. Ilgesnes nei 100 simboliu
    """

    filter_long_fragments(seq_frames)
    print('\n\n')
    filter_long_fragments(reverse_seq_frames, reverse_seq = True)
    print('\n\n')

    """
    4. Parasykite funkcijas, kurios ivertintu kodonu ir dikodonu daznius
    (visi imanomi kodonai/dikodonai ir jų atitinkamas daznis  - gali buti nemazai nuliu, jei ju sekoje nerasite).
    """

    get_frames_codon_freq(seq_frames)
    get_frames_codon_freq(reverse_seq_frames)

    get_frames_decodon_freq(seq_frames)
    get_frames_decodon_freq(reverse_seq_frames)

    get_total_codon_freq(seq_codons, show_table = True)
    get_total_decodon_freq(seq_decodons, show_table = True)

if __name__ == "__main__":
    main()