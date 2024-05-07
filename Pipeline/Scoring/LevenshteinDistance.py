from Levenshtein import distance

from Pipeline.Scoring.AScore import AScore


class LevenshteinDistance(AScore):
    def __init__(self):
        pass

    def getScore(self, predicted:str, actual:str) ->float:
        """Calculate the Levenshtein distance between two sequences.
        Parameters:
        :param predicted: First (predicted) sequence.
        :param actual: Second (actual) sequence.
        :return: Levenshtein distance between the two sequences.
        """
        return distance(predicted, actual)

if __name__ == "__main__":
    levenshtein = LevenshteinDistance()
    print(levenshtein.getScore("ITHQEVDSR", "LTHQEVDSR"))
    print(levenshtein.getScore(actual="VAMAMGSHIR", predicted="GSHP"))
    print(levenshtein.getScore(predicted="GSHL", actual="VAMAMGSHIR"))