import pandas as pd

from Pipeline.AlgorithmResultParsers.DeepNovoParser import DeepNovoParser
from Pipeline.AlgorithmResultParsers.NovorParser import NovorParser
from Pipeline.Scoring.SequenceIdentity import SequenceIdentity
from Pipeline.Scoring.SequenceSimilarity import SequenceSimilarity


def cleanPeptideModifications(peptide: str) -> str:
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
    peptide = entry['output_seq']
    # peptide = cleanPeptideModifications(entry['peptide'])
    return [peptide, entry['Sequence'], similarity.getScore(peptide, entry['Sequence']),
            identiy.getScore(peptide, entry['Sequence'])]


if __name__ == "__main__":
    """
    novor = NovorParser(
        "/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/BD7_Thermo_Pool52_HCD/Novor/Run_2/01640c_BD7-Thermo_SRM_Pool_52_01_01-3xHCD-1h-R2.mzml.novor.csv",
        20)
    expected = pd.read_csv(
        "/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/BD7_Thermo_Pool52_HCD/Thermo_SRM_Pool_52_01_01_3xHCD-1h-R2-tryptic/msmsScans.txt",
        sep='\t')
    expected_seq = pd.DataFrame(expected['Sequence'], )
    merged_df = pd.merge(novor.getParsedData(), expected_seq, left_on='id', right_index=True, how='inner')
    # print(merged_df)
    res = merged_df.query('Sequence != \' \'').apply(lambda x: calculateScores(x), axis=1)
    analysis = pd.DataFrame(list(res), columns=['Peptide', 'Sequence', 'Similarity', 'Identity'])
    analysis.to_csv(
        "/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/Result/Analysis_Novor.tsv",
        sep='\t')
"""
    deepnovo_df = DeepNovoParser(
        "/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/Result/Analysis_DeepNovo.tsv")
    expected = pd.read_csv(
        "/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/BD7_Thermo_Pool52_HCD/Thermo_SRM_Pool_52_01_01_3xHCD-1h-R2-tryptic/msmsScans.txt",
        sep='\t')
    expected_seq = pd.DataFrame({'Scan number': expected['Scan number'], 'Sequence': expected['Sequence']})
    merged_df = pd.merge(deepnovo_df.getParsedData(), expected_seq, left_on='scan', right_on='Scan number', how='inner')
    print(merged_df)
    res = merged_df.query('Sequence != \' \'').apply(lambda x: calculateScores(x), axis=1)
    analysis = pd.DataFrame(list(res), columns=['Peptide', 'Sequence', 'Similarity', 'Identity'])
    print(analysis)
    analysis.to_csv(
        "/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/Result/Analysis_DeepNovo_Result.tsv",
        sep='\t')
