from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

from Pipeline.Scoring.AScore import AScore

class SequenceSimilarity(AScore):
    def __init__(self, substitution_matrix:str = 'BLOSUM62', alignment_mode:str = 'global', open_gap_score:int = -2, extend_gap_score:int = -2):
        self.substitution_matrix = substitution_matrix
        self.aligner = PairwiseAligner()
        self.aligner.mode = alignment_mode
        self.aligner.open_gap_score = open_gap_score
        self.aligner.extend_gap_score = extend_gap_score
        self.aligner.substitution_matrix = substitution_matrices.load(self.substitution_matrix)
    def getScore(self, predicted:str, actual:str)->float:
        """Calculate the percent similarity between two sequences.
        Parameters:
        :param predicted: First (predicted) sequence.
        :param actual: Second (actual) sequence.
        :return: Percent similarity between the two sequences.
        """

        # align the sequences
        alignments = self.aligner.align(predicted, actual)
        # get the first (best) alignment
        if len(alignments)==0:
            return 0.0
        alignment = alignments[0]
        # calculate the number of positive substitution scores
        positive_substitution_scores = sum(
            1 for a, b in zip(alignment[0,:], alignment[1,:]) if self.aligner.substitution_matrix.get((a, b),self.aligner.open_gap_score) >0)
        # precent similarity is the number of positive substitution scores divided by the length of the alignment
        length = len(alignment[0,:]) if self.aligner.mode == 'global' else len(predicted)
        similarity = positive_substitution_scores/length

        return similarity


if __name__ == "__main__":
    similarity = SequenceSimilarity()
    print(similarity.getScore("ITHQGEVDSR","LTHQEPVDSR"))
    similarity = SequenceSimilarity(alignment_mode='local')
    print(similarity.getScore(actual="AFPSPQTLLEDPLR", predicted="PEKD"))
    print(similarity.getScore(predicted="PYGC", actual="SHTGEKPYGCNECGK"))