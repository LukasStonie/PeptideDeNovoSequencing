from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

from Pipeline.Scoring.AScore import AScore

class SequenceIdentity(AScore):
    def __init__(self, substitution_matrix:str = 'BLOSUM62', alignment_mode:str = 'global', open_gap_score:int = -2, extend_gap_score:int = -2):
        self.substitution_matrix = substitution_matrix
        self.aligner = PairwiseAligner()
        self.aligner.mode = alignment_mode
        self.aligner.open_gap_score = open_gap_score
        self.aligner.extend_gap_score = extend_gap_score
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
        print(alignment)
        # calculate the number of identical positions
        identical_positions = sum(a==b for a,b in zip(alignment[0,:], alignment[1,:]))

        # precent identity is the number of identical positions divided by the length of the alignment
        length =  len(alignment[0,:]) if self.aligner.mode == 'global' else len(predicted)
        identity = identical_positions/length

        return identity

if __name__ == "__main__":
    identity = SequenceIdentity()
    print(identity.getScore(predicted="FELATVTEK", actual="FQIATVTEK"))
    identity = SequenceIdentity(alignment_mode='local', open_gap_score=-10, extend_gap_score=-10)
    print(identity.getScore(actual="LEESLATTETFK", predicted="EEST"))
    print(identity.getScore(predicted="GSHT", actual="VAMAMGSHPR"))