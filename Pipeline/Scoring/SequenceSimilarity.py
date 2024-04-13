from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
class SequenceSimilarity():
    def __init__(self, substitution_matrix:str = 'BLOSUM62'):
        self.substitution_matrix = substitution_matrix
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.open_gap_score = -2
        self.aligner.extend_gap_score = -2
        self.aligner.substitution_matrix = substitution_matrices.load(self.substitution_matrix)
    def getScore(self, sequence1:str, sequence2:str)->float:
        """Calculate the percent similarity between two sequences.
        Parameters:
        :param sequence1: First (predicted) sequence.
        :param sequence2: Second (actual) sequence.
        :return: Percent similarity between the two sequences.
        """

        # align the sequences
        alignments = self.aligner.align(sequence1, sequence2)
        # get the first (best) alignment
        alignment = alignments[0]
        # calculate the number of positive substitution scores
        positive_substitution_scores = sum(
            1 for a, b in zip(alignment[0,:], alignment[1,:]) if self.aligner.substitution_matrix.get((a, b),self.aligner.open_gap_score) >0)
        # precent similarity is the number of positive substitution scores divided by the length of the alignment
        length = len(alignment[0,:])
        similarity = positive_substitution_scores/length

        return similarity


if __name__ == "__main__":
    similarity = SequenceSimilarity()
    print(similarity.getScore("ITHQGEVDSR","LTHQEPVDSR"))