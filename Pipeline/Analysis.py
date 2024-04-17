import pandas as pd

from Pipeline.AlgorithmResultParsers.DeepNovoParser import DeepNovoParser
from Pipeline.AlgorithmResultParsers.NovorParser import NovorParser
from Pipeline.Scoring.SequenceIdentity import SequenceIdentity
from Pipeline.Scoring.SequenceSimilarity import SequenceSimilarity


def calculateScores(entry)->list:
    """Calculate the similarity and identity of two sequences.
    Parameters:
    :param entry: DataFrame entry containing the predicted and actual sequence.
    :return: List containing the predicted sequence, the actual sequence, the similarity score, and the identity score.
        """
    identiy = SequenceIdentity()
    similarity = SequenceSimilarity()
    return [entry['Predicted'], entry['Actual'], similarity.getScore(entry['Predicted'], entry['Actual']),
            identiy.getScore(entry['Predicted'], entry['Actual'])]


if __name__ == "__main__":
    # read expected sequences
    expected = pd.read_csv(
        "/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/BD7_Thermo_Pool52_HCD/Thermo_SRM_Pool_52_01_01_3xHCD-1h-R2-tryptic/msmsScans.txt",
        sep='\t')
    actualSequence_df = pd.DataFrame({'Scan': expected['Scan number'], 'Actual': expected['Sequence']})

    # read Novor result
    novor_df = NovorParser(
        "/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/BD7_Thermo_Pool52_HCD/Novor/Run_2/01640c_BD7-Thermo_SRM_Pool_52_01_01-3xHCD-1h-R2.mzml.novor.csv",
        20).parse()
    # merge expected and actual sequences for novor
    merged_df = pd.merge(novor_df, actualSequence_df, left_on='ID', right_index=True, how='inner')
    # score the sequences and write the result to a file
    novor_result = merged_df.query('Actual != \' \'').apply(lambda x: calculateScores(x), axis=1)
    novor_result_df = pd.DataFrame(list(novor_result), columns=['Predicted', 'Actual', 'Similarity', 'Identity'])
    novor_result_df.to_csv(
        "/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/Result/Analysis_Novor_1.tsv",
        sep='\t')

    # read DeepNovo result
    deepnovo_df = DeepNovoParser(
        "/Data/BD7_Thermo_Pool52_HCD/Result_run_1/Result_run_1.tsv").parse()
    # merge expected and actual sequences for deepnovo
    merged_df = pd.merge(deepnovo_df, actualSequence_df, left_on='Scan', right_on='Scan', how='inner')
    # score the sequences and write the result to a file
    deepnovo_result = merged_df.query('Actual != \' \'').apply(lambda x: calculateScores(x), axis=1)
    deepnovo_result_df = pd.DataFrame(list(deepnovo_result), columns=['Peptide', 'Sequence', 'Similarity', 'Identity'])
    deepnovo_result_df.to_csv(
        "/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/Result/Analysis_DeepNovo_1.tsv",
        sep='\t')
