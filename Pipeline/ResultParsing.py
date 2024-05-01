import pandas as pd

from Pipeline.AlgorithmResultParsers.DeepNovoParser import DeepNovoParser
from Pipeline.AlgorithmResultParsers.DirecTagParser import DirecTagParser
from Pipeline.AlgorithmResultParsers.NovorParser import NovorParser
from Pipeline.AlgorithmResultParsers.PEAKSParser import PEAKSParser


def readPeaksResult(file: str, actual_sequences: pd.DataFrame) -> pd.DataFrame:
    """Reads the PEAKS result file.
    Parameters:
    :param file: Path to the PEAKS result file.
    :return: DataFrame containing the PEAKS result merged with the actual sequences.
    """
    # read PEAKS result
    peaks_df = PEAKSParser(file).parse()
    # merge expected and actual sequences for peaks
    merged_df = pd.merge(peaks_df, actual_sequences, left_on='Scan', right_on='Scan', how='inner')
    # score the sequences and write the result to a file
    peaks_result = merged_df.query('Actual != \' \'')
    return peaks_result


def readDirecTagResult(file: str, actual_sequences: pd.DataFrame) -> pd.DataFrame:
    """Reads the DirecTag result file.
    Parameters:
    :param file: Path to the DirecTag result file.
    :return: DataFrame containing the DirecTag result merged with the actual sequences.
    """
    # read DirecTag result
    direcTag_df = DirecTagParser(file, 25).parse()
    # merge expected and actual sequences for direcTag
    merged_df = pd.merge(direcTag_df, actual_sequences, left_on='ID', right_index=True, how='inner')
    # score the sequences and write the result to a file
    direcTag_result = merged_df.query('Actual != \' \'')
    return direcTag_result


def readNovorResult(file: str, actual_sequences: pd.DataFrame) -> pd.DataFrame:
    # read Novor result
    novor_df = NovorParser(file, 20).parse()
    # merge expected and actual sequences for novor
    merged_df = pd.merge(novor_df, actual_sequences, left_on='ID', right_index=True, how='inner')
    # score the sequences and write the result to a file
    novor_result = merged_df.query('Actual != \' \'')
    return novor_result


def readDeepNovoResult(file: str, actual_sequences: pd.DataFrame) -> pd.DataFrame:
    # read DeepNovo result
    deepnovo_df = DeepNovoParser(file).parse()
    # merge expected and actual sequences for deepnovo
    merged_df = pd.merge(deepnovo_df, actual_sequences, left_on='Scan', right_on='Scan', how='inner')
    # score the sequences and write the result to a file
    deepnovo_result = merged_df.query('Actual != \' \'')
    return deepnovo_result


if __name__ == "__main__":
    # read expected sequences for pool 49

    pools = [
        ('Pool_49', 'BA7'),
        ('Pool_52', 'BD7'),
        ('Pool_60', 'BD8')]

    for p in pools:
        expected_pool_49 = pd.read_csv(
            f"../Data/Datasets/{p[0]}/Thermo_SRM_{p[0]}_01_01_3xHCD-1h-R2-tryptic/msmsScans.txt",
            sep='\t')
        actualSequence_pool_df = pd.DataFrame(
            {'Scan': expected_pool_49['Scan number'], 'Actual': expected_pool_49['Sequence']})

        # read peaks result
        peaks_result = readPeaksResult(f"../Data/AlgorithmResults/{p[0]}/PEAKS/Sample 1.denovo.csv",
                                       actualSequence_pool_df)
        # write the result to a file
        peaks_result.to_csv(f'../Data/ParsingResults/{p[0]}/peaks_results.tsv', sep='\t', index=None)

        # read direcTag result
        direcTag_result = readDirecTagResult(
            f"../Data/AlgorithmResults/{p[0]}/DirecTag/Run_1/01640c_{p[1]}-Thermo_SRM_{p[0]}_01_01-3xHCD-1h-R2.tags",
            actualSequence_pool_df)
        # write the result to a file
        direcTag_result.to_csv(f'../Data/ParsingResults/{p[0]}/direcTag_results.tsv', sep='\t', index=None)

        # read novor result
        novor_result = readNovorResult(
            f"../Data/AlgorithmResults/{p[0]}/Novor/Run_1/01640c_{p[1]}-Thermo_SRM_{p[0]}_01_01-3xHCD-1h-R2.novor.csv",
            actualSequence_pool_df)
        # write the result to a file
        novor_result.to_csv(f'../Data/ParsingResults/{p[0]}/novor_results.tsv', sep='\t', index=None)

        # read deepnovo result
        deepnovo_result = readDeepNovoResult(f"../Data/AlgorithmResults/{p[0]}/DeepNovo/decode_output.tab",
                                             actualSequence_pool_df)
        # write the result to a file
        deepnovo_result.to_csv(f'../Data/ParsingResults/{p[0]}/deepnovo_results.tsv', sep='\t', index=None)
