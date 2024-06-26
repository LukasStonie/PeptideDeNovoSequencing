from Bio.Align import substitution_matrices, PairwiseAligner

from Pipeline.Scoring.AScore import AScore


class AlignmentScore(AScore):
    def __init__(self, substitution_matrix:str = 'BLOSUM62', alignment_mode:str = 'global', open_gap_score:int = -2, extend_gap_score:int = -2):
        self.substitution_matrix = substitution_matrix
        self.aligner = PairwiseAligner()
        self.aligner.mode = alignment_mode
        self.aligner.open_gap_score = open_gap_score
        self.aligner.extend_gap_score = extend_gap_score
        self.aligner.substitution_matrix = substitution_matrices.load(self.substitution_matrix)

    def getScore(self, predicted: str, actual: str) -> float:
        """Perform a global alignment between two sequences and calculate the alignment score.
        Parameters:
        :param predicted: First sequence.
        :param actual: Second sequence.
        :return: Alignment score between the two sequences.
            """

        # align the sequences
        alignments = self.aligner.align(predicted, actual)
        if len(alignments)==0:
            return 0.0
        # get the first (best) alignment
        alignment = alignments[0]
        # return the alignment score
        return alignment.score


if __name__ == "__main__":
    alignment = AlignmentScore()
    print(alignment.getScore("ITHQGEVDSR", "LTHQEVDSR"))
    alignment = AlignmentScore(alignment_mode='local')
    print(alignment.getScore(actual="VAMAMGSHPR", predicted="QWTY"))