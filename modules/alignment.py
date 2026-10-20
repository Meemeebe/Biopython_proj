
from Bio.Align import PairwiseAligner


def align_sequences(sequence1, sequence2):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    return aligner.align(
        str(sequence1).upper(),
        str(sequence2).upper()
    )[0]
