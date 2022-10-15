from utils.process_virus import get_virus, get_total_codon_freq, get_total_decodon_freq
from utils.utils import get_distance_matrix, write_to_phylip

"""
5. Palyginkite kodonu bei dikodonu daznius tarp visu seku
"""

bacterial1 = get_virus('bacterial1.fasta')
bacterial2 = get_virus('bacterial2.fasta')
bacterial3 = get_virus('bacterial3.fasta')
bacterial4 = get_virus('bacterial4.fasta')

mamalian1 = get_virus('mamalian1.fasta')
mamalian3 = get_virus('mamalian2.fasta')
mamalian2 = get_virus('mamalian3.fasta')
mamalian4 = get_virus('mamalian4.fasta')

virus_list = [bacterial1, bacterial2, bacterial3, bacterial4, mamalian1, mamalian2, mamalian3, mamalian4]

virus_codon_totals = {virus.id : virus.codon_freq for virus in virus_list}
get_total_codon_freq(virus_codon_totals)
codon_dist_matrix = get_distance_matrix(virus_codon_totals)
print(f'Atstumu matrica remiantis kodonu dazniu lentele:\n\n{codon_dist_matrix}\n')


virus_decodon_totals = {virus.id : virus.decodon_freq for virus in virus_list}
get_total_decodon_freq(virus_decodon_totals)
decodon_dist_matrix = get_distance_matrix(virus_decodon_totals)
print(f'Atstumu matrica remiantis dekodonu dazniu lentele:\n\n{decodon_dist_matrix}')

write_to_phylip('codon_dist_matrix.phy', codon_dist_matrix)
write_to_phylip('decodon_dist_matrix.phy', decodon_dist_matrix)
