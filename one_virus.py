from utils.utils import *
from virus_structure import *
from virus_parts import *
import sys

def main():
    """
    1. Pateiktoje sekoje fasta formatu surastu visas start ir stop kodonų poras,
    tarp kurių nebutu stop kodono (ir tiesioginei sekai ir jos reverse komplementui). 
    2. Kiekvienam stop kodonui parinkti toliausiai nuo jo esanti start kodoną (su salyga, kad tarp ju nera kito stop kodono)
    """

    filename = sys.argv[1]
    virus = read_file(filename)
    # virus_sequence = Seq('AAATTAG')

    virus_sequence = virus.seq
    reverse_virus_sequence = reverse_sequence(virus_sequence)

            # j = name[key]['start_stop_p']
            # print(f'Seka: {straight}; skaitymo remelis: {key}. Start ir stop kodonų poros: {j}\n')

    part1(virus_sequence, seq_frames)
    print()
    part1(reverse_virus_sequence, reverse_seq_frames)


    """
    3. Atfiltruokite visus fragmentus ("tai butu baltymų koduojancios sekos"), kurie trumpesni nei 100 fragmentų. Ilgesnes nei 100 simboliu
    """

    def filter_long_fragments(frames):
        for key in frames:
            long_frag = get_long_fragments(frames[key]['seq'], frames[key]['start_stop_p'])
            # print(f'Seka: {straight}; skaitymo remelis: {key}. Fragmentai ilgesni nei 100 simboliu: {long_frag}\n') 


    filter_long_fragments(seq_frames)
    print('\n')
    filter_long_fragments(reverse_seq_frames)
    print('\n')


    """
    4. Parasykite funkcijas, kurios ivertintu kodonu ir dikodonu daznius
    (visi imanomi kodonai/dikodonai ir jų atitinkamas daznis  - gali buti nemazai nuliu, jei ju sekoje nerasite).
    """

    part2(seq_frames)
    part2(reverse_seq_frames)

    part3(seq_frames)
    part3(reverse_seq_frames)

    part4(seq_codons, show_table = True)
    part5(seq_decodons, show_table = True)

if __name__ == "__main__":
    main()