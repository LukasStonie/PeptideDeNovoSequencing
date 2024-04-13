import pandas as pd

from Pipeline.AlgorithmResultParsers.NovorParser import NovorParser
from Pipeline.Scoring.SequenceIdentity import SequenceIdentity
from Pipeline.Scoring.SequenceSimilarity import SequenceSimilarity

def cleanPeptideModifications(peptide:str)->str:
    """Remove all modifications from a peptide sequence.
    Parameters:
    :param peptide: Peptide sequence with modifications.
    :return: Peptide sequence without modifications.
    """
    return ''.join([aa for aa in peptide if aa.isalpha()])
def calculateScores(entry):
    if entry['Sequence'] == ' ':
        return
    identiy = SequenceIdentity()
    similarity = SequenceSimilarity()
    peptide = cleanPeptideModifications(entry['peptide'])
    print(peptide, entry['Sequence'], similarity.getScore(peptide, entry['Sequence']), identiy.getScore(peptide, entry['Sequence']) )

if __name__ == "__main__":
    novor = NovorParser("/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/BD7_Thermo_Pool52_HCD/Novor/Run_2/01640c_BD7-Thermo_SRM_Pool_52_01_01-3xHCD-1h-R2.mzml.novor.csv", 20)
    expected = pd.read_csv("/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/BD7_Thermo_Pool52_HCD/Thermo_SRM_Pool_52_01_01_3xHCD-1h-R2-tryptic/msmsScans.txt", sep='\t')
    expected_seq = pd.DataFrame(expected['Sequence'])
    merged_df = pd.merge(novor.getParsedData(), expected_seq, left_on='id', right_index=True, how='inner')
    # print(merged_df)
    merged_df.head(100).apply(lambda x: calculateScores(x), axis=1)