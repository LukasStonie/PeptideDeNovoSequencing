from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
class SequenceIdentity():
    def __init__(self, substitution_matrix:str = 'PAM250'):
        self.substitution_matrix = substitution_matrix
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.open_gap_score = -10
        self.aligner.extend_gap_score = -0.5
        self.aligner.substitution_matrix = substitution_matrices.load(self.substitution_matrix)
    def getScore(self, sequence1:str, sequence2:str)->float:
        """Calculate the percent identity between two sequences.
        Parameters:
        :param sequence1: First (predicted) sequence.
        :param sequence2: Second (actual) sequence.
        :return: Percent identity between the two sequences.
        """
        # align the sequences
        alignments = self.aligner.align(sequence1, sequence2)

        # get the first (best) alignment
        alignment = alignments[0]

        # calculate the number of identical positions
        identical_positions = sum(a==b for a,b in zip(alignment[0,:], alignment[1,:]))

        # precent identity is the number of identical positions divided by the length of the alignment
        length = len(alignment[0,:])
        identity = identical_positions/length

        return identity

if __name__ == "__main__":
    similarity = SequenceIdentity()
    print(similarity.getScore("ITHQGEVDSR","LTHQEVDSR"))