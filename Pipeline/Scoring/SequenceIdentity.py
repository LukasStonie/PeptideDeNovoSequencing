from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

from Pipeline.Scoring.AScore import AScore

class SequenceIdentity(AScore):
    def __init__(self, substitution_matrix:str = 'PAM250', alignment_mode:str = 'global'):
        self.substitution_matrix = substitution_matrix
        self.aligner = PairwiseAligner()
        self.aligner.mode = alignment_mode
        self.aligner.open_gap_score = -10
        self.aligner.extend_gap_score = -0.5
        self.aligner.substitution_matrix = substitution_matrices.load(self.substitution_matrix)
    def getScore(self, predicted:str, actual:str)->float:
        """Calculate the percent identity between two sequences.
        Parameters:
        :param predicted: First (predicted) sequence.
        :param actual: Second (actual) sequence.
        :return: Percent identity between the two sequences.
        """
        # align the sequences
        alignments = self.aligner.align(predicted, actual)

        if len(alignments)==0:
            return 0.0
        # get the first (best) alignment
        alignment = alignments[0]

        # calculate the number of identical positions
        identical_positions = sum(a==b for a,b in zip(alignment[0,:], alignment[1,:]))

        # precent identity is the number of identical positions divided by the length of the alignment
        length =  len(alignment[0,:]) if self.aligner.mode == 'global' else len(predicted)
        identity = identical_positions/length

        return identity

if __name__ == "__main__":
    identity = SequenceIdentity()
    print(identity.getScore("ITHQGEVDSR", "LTHQEVDSR"))
    identity = SequenceIdentity(alignment_mode='local')
    print(identity.getScore(actual="VAMAMGSHPR", predicted="GSHP"))
    print(identity.getScore(predicted="GSHT", actual="VAMAMGSHPR"))