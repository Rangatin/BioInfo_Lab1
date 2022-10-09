from virus_structure import *
from utils import *


def part1(virus_sequence, name):
    i = 0
    for key in name:
        name[key]['seq'] = virus_sequence[i:]
        name[key]['codons'] = get_codons(name[key]['seq'])
        name[key]['start_stop_p'] = get_start_stop_pairs(name[key]['codons'])
        i += 1


def part2(virus_frames):
    for key in virus_frames:
        codon_freq = Counter(virus_frames[key]['codons'])
        l = seq_codons[key]
        for el in l:
            l[el] = codon_freq[el]


def part3(virus_frames):
    for key in virus_frames:
        virus_frames[key]['decodon_freq'] = Counter(get_decodons(virus_frames[key]['seq']))
        l = seq_decodons[key]
        for el in l:
            l[el] = virus_frames[key]['decodon_freq'][el]


def part4(seq_codons, show_table = False):
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


def part5(seq_decodons, show_table = False):
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
    virus_seq = self.seq
    reverse_virus_sequence = reverse_sequence(virus_seq)
    virus_frames = seq_frames
    reverse_virus_frames = reverse_seq_frames

    part1(virus_seq, virus_frames)
    part1(reverse_virus_sequence, reverse_virus_frames)

    part2(virus_frames)
    part2(reverse_virus_frames)

    part3(virus_frames)
    part3(virus_frames)

    total_codon_freq = part4(seq_codons)
    total_decodon_freq = part5(seq_decodons)
    return (total_codon_freq, total_decodon_freq)