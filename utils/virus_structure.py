from utils.utils import get_all_possible_codons, get_all_possible_decodons

def get_blueprint():
    frame_blueprint = {
        'seq' : [],
        'codons': [],
        'start_stop_p' : [],
        # 'long_frag' : [],
        'codon_freq' : get_all_possible_codons(),
        'decodon_freq' : get_all_possible_decodons(),
    }
    return frame_blueprint


seq_frames = {
    'frame1' : get_blueprint(),
    'frame2' : get_blueprint(),
    'frame3' : get_blueprint(),
}

reverse_seq_frames = {
    'reverse_frame1' : get_blueprint(),
    'reverse_frame2' : get_blueprint(),
    'reverse_frame3' : get_blueprint(),
}



seq_codons = {
    'frame1' : get_all_possible_codons(),
    'frame2' : get_all_possible_codons(),
    'frame3' : get_all_possible_codons(),
    'reverse_frame1' : get_all_possible_codons(),
    'reverse_frame2' : get_all_possible_codons(),
    'reverse_frame3' : get_all_possible_codons(),
}

seq_decodons = {
    'frame1' : get_all_possible_decodons(),
    'frame2' : get_all_possible_decodons(),
    'frame3' : get_all_possible_decodons(),
    'reverse_frame1' : get_all_possible_decodons(),
    'reverse_frame2' : get_all_possible_decodons(),
    'reverse_frame3' : get_all_possible_decodons(),
}

