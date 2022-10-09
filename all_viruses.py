import pandas as pd
from utils.utils import read_file
from virus_parts import get_virus_info, part5, part4
from scipy.spatial.distance import pdist, squareform

"""
5. Palyginkite kodonu bei dikodonu daznius tarp visu seku
"""
def get_virus(filename):
    virus = read_file(filename)
    (virus.codon_freq, virus.decodon_freq) = get_virus_info(virus)
    return virus

bacterial1 = get_virus('bacterial1.fasta')
bacterial2 = get_virus('bacterial2.fasta')
bacterial3 = get_virus('bacterial3.fasta')
bacterial4 = get_virus('bacterial4.fasta')

mamalian1 = get_virus('mamalian1.fasta')
mamalian3 = get_virus('mamalian2.fasta')
mamalian2 = get_virus('mamalian3.fasta')
mamalian4 = get_virus('mamalian4.fasta')

oa = [bacterial1, bacterial2, bacterial3, bacterial4, mamalian1, mamalian2, mamalian3, mamalian4]

def get_distance_matrix(dict):
    virus_df = pd.DataFrame(dict).T

    pairwise = pd.DataFrame(
        squareform(pdist(virus_df)),
        columns = virus_df.index,
        index = virus_df.index
    )

    return pairwise

virus_codon_totals = {virus.id : virus.codon_freq for virus in oa}
part4(virus_codon_totals)
codon_dist_matrix = get_distance_matrix(virus_codon_totals)

print(f'Atstumu matrica remiantis kodonu dazniu lentele:\n\n{codon_dist_matrix}\n')

virus_decodon_totals = {virus.id : virus.decodon_freq for virus in oa}
part5(virus_decodon_totals)
decodon_dist_matrix = get_distance_matrix(virus_decodon_totals, key='decodon')

print(f'Atstumu matrica remiantis dekodonu dazniu lentele:\n\n{decodon_dist_matrix}')

with open('codon_dist_matrix.txt', 'w') as f:
    dfAsString = codon_dist_matrix.to_string(header=False)
    f.write(dfAsString)