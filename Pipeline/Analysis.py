import pandas as pd
from joblib import Parallel, delayed


from Pipeline.AlgorithmResultParsers.DeepNovoParser import DeepNovoParser
from Pipeline.AlgorithmResultParsers.DirecTagParser import DirecTagParser
from Pipeline.AlgorithmResultParsers.NovorParser import NovorParser
from Pipeline.Scoring.AlignmentScore import AlignmentScore
from Pipeline.Scoring.LevenshteinDistance import LevenshteinDistance
from Pipeline.Scoring.NormalizedAlignmentScore import NormalizedAlignmentScore
from Pipeline.Scoring.SequenceIdentity import SequenceIdentity
from Pipeline.Scoring.SequenceSimilarity import SequenceSimilarity

RESULT_CLOUMNS = ['ID','Predicted', 'Actual', 'Similarity', 'Identity', 'Local Alignment', 'Global Alignment','Normalized Local Alignment', 'Normalized Global Alignment', 'Levenshtein']

def calculateScores(entry, alignment_mode:str = 'global') -> list:
    """Calculate the similarity and identity of two sequences.
    Parameters:
    :param entry: DataFrame entry containing the predicted and actual sequence.
    :return: List containing the predicted sequence, the actual sequence, the similarity score, and the identity score.
        """
    identiy = SequenceIdentity(alignment_mode=alignment_mode)
    similarity = SequenceSimilarity(alignment_mode=alignment_mode)
    localAlignmentScore = AlignmentScore(alignment_mode='local')
    globalAlignmentScore = AlignmentScore(alignment_mode='global')
    normalizedLocalAlignmentScore = NormalizedAlignmentScore(alignment_mode='local')
    normalizedGlobalAlignmentScore = NormalizedAlignmentScore(alignment_mode='global')
    levenshtein = LevenshteinDistance()
    return [
        entry['ID'] if 'ID' in entry else entry['Scan'],
        entry['Predicted'], entry['Actual'],
        similarity.getScore(predicted=entry['Predicted'], actual=entry['Actual']),
        identiy.getScore(predicted=entry['Predicted'],actual= entry['Actual']),
        localAlignmentScore.getScore(predicted=entry['Predicted'],actual= entry['Actual']),
        globalAlignmentScore.getScore(predicted=entry['Predicted'],actual= entry['Actual']),
        normalizedLocalAlignmentScore.getScore(predicted=entry['Predicted'], actual=entry['Actual']),
        normalizedGlobalAlignmentScore.getScore(predicted=entry['Predicted'], actual=entry['Actual']),
        levenshtein.getScore(predicted=entry['Predicted'],actual= entry['Actual'])]


def processNovor(file: str, actualSequence_df: pd.DataFrame):
    # read Novor result
    novor_df = NovorParser(file, 20).parse()
    # merge expected and actual sequences for novor
    merged_df = pd.merge(novor_df, actualSequence_df, left_on='ID', right_index=True, how='inner')
    # score the sequences and write the result to a file
    novor_result = merged_df.query('Actual != \' \'').apply(lambda x: calculateScores(x), axis=1)
    return pd.DataFrame(list(novor_result), columns=RESULT_CLOUMNS)

def processDeepNovo(file: str, actualSequence_df: pd.DataFrame):
    # read DeepNovo result
    deepnovo_df = DeepNovoParser(file).parse()
    # merge expected and actual sequences for deepnovo
    merged_df = pd.merge(deepnovo_df, actualSequence_df, left_on='Scan', right_on='Scan', how='inner')
    # score the sequences and write the result to a file
    deepnovo_result = merged_df.query('Actual != \' \'').apply(lambda x: calculateScores(x), axis=1)
    return pd.DataFrame(list(deepnovo_result), columns=RESULT_CLOUMNS)


def calculateChunkIndex(len_df, num_cores):
    """Calculate the start and end indices for each chunk.
    Parameters:
    :param len_df: Length of the DataFrame.
    :param num_cores: Number of cores to use for parallel processing.
    :return: List of tuples containing the start and end indices for each chunk.
    """
    # Calculate the chunk size
    chunk_size = len_df // num_cores
    # Initialize a list to store start and end indices for each chunk
    chunk_indices = []
    # Create chunks and calculate start and end indices
    for i in range(num_cores):
        start_index = i * chunk_size
        end_index = (i + 1) * chunk_size if i < num_cores - 1 else len_df
        chunk_indices.append((start_index, end_index))

    return chunk_indices
def processDirecTag(subset_df, alignment_mode:str = 'global'):
    output = subset_df.apply(lambda x: calculateScores(x, alignment_mode), axis=1)
    output = pd.DataFrame(list(output), columns=RESULT_CLOUMNS)
    return output
def processDirecTag_parallel(file:str, actualSequence_df:pd.DataFrame):
    # read DirecTag result
    direcTag_df = DirecTagParser(file, 25).parse()
    # merge expected and actual sequences for direcTag
    merged_df = pd.merge(direcTag_df, actualSequence_df, left_on='ID', right_index=True, how='inner')
    # score the sequences and write the result to a file
    direcTag_result = merged_df.query('Actual != \' \'')

    # define number of cpus and calculate chunk indices for parallel processing
    cpus = 6  # os.cpu_count()
    chunk_indices = calculateChunkIndex(len(direcTag_result), cpus)

    # define parallel processing object with number of cpus
    parallel_direcTag = Parallel(n_jobs=cpus)
    # process the direcTag result in parallel
    output = parallel_direcTag(delayed(processDirecTag)(direcTag_result.iloc[chunk[0]:chunk[1], :], 'local') for chunk in chunk_indices)
    # concatenate the results
    output = pd.concat(output, axis=0)
    return output

if __name__ == "__main__":
    # read expected sequences
    expected = pd.read_csv(
        "/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/BD7_Thermo_Pool52_HCD/Thermo_SRM_Pool_52_01_01_3xHCD-1h-R2-tryptic/msmsScans.txt",
        sep='\t')
    actualSequence_df = pd.DataFrame({'Scan': expected['Scan number'], 'Actual': expected['Sequence']})


    novor_df = processNovor(
        "/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/BD7_Thermo_Pool52_HCD/Novor/Run_2/01640c_BD7-Thermo_SRM_Pool_52_01_01-3xHCD-1h-R2.mzml.novor.csv",
        actualSequence_df)
    novor_df.to_csv(
        "/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/Result/Analysis_Novor_1.tsv",
        sep='\t', index=False)

    deepnovo_df = processDeepNovo("/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/BD7_Thermo_Pool52_HCD/DeepNovo/Result_run_1.tsv", actualSequence_df)
    deepnovo_df.to_csv(
        "/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/Result/Analysis_DeepNovo_1.tsv",
        sep='\t', index=False)

    direcTag_df = processDirecTag_parallel("/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/BD7_Thermo_Pool52_HCD/DirecTag/Run_1/01640c_BD7-Thermo_SRM_Pool_52_01_01-3xHCD-1h-R2.tags",
                                           actualSequence_df)
    print(direcTag_df)
    direcTag_df.to_csv("/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/Result/Analysis_DirecTag_1.tsv",sep='\t', index=False)