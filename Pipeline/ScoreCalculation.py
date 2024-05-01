import pandas as pd
from joblib import Parallel, delayed

from Scoring.AlignmentScore import AlignmentScore
from Scoring.LevenshteinDistance import LevenshteinDistance
from Scoring.NormalizedAlignmentScore import NormalizedAlignmentScore
from Scoring.SequenceIdentity import SequenceIdentity
from Scoring.SequenceSimilarity import SequenceSimilarity

import warnings

warnings.filterwarnings("ignore")

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
def calculateScoresOfChunk(subset_df, alignment_mode:str = 'global'):
    output = subset_df.apply(lambda x: calculateScores(x, alignment_mode), axis=1)
    output = pd.DataFrame(list(output), columns=RESULT_CLOUMNS)
    return output

def processParsed(file:str):
    # read parsed result
    parsed_df = pd.read_csv(file, sep='\t', index_col=None, header=0)

    # define number of cpus and calculate chunk indices for parallel processing
    cpus = 6  # os.cpu_count()
    chunk_indices = calculateChunkIndex(len(parsed_df), cpus)

    # define parallel processing object with number of cpus
    parallel_direcTag = Parallel(n_jobs=cpus)
    # process the direcTag result in parallel
    output = parallel_direcTag(delayed(calculateScoresOfChunk)(parsed_df.iloc[chunk[0]:chunk[1], :], 'local') for chunk in chunk_indices)
    # concatenate the results
    output = pd.concat(output, axis=0)
    return output

if __name__ == "__main__":
    pools = ['Pool_49','Pool_52','Pool_60']

    for p in pools:
        print("Scoring ", p, " - PEAKS")
        # score peaks
        peaks_scored_df = processParsed(f"../Data/ParsingResults/{p}/peaks_results.tsv")
        # save peaks scores
        peaks_scored_df.to_csv(f"../Data/ScoringResults/{p}/peaks_scored.tsv", sep='\t', index=None)

        print("Scoring ", p, " - Novor")
        # score novor
        novor_scored_df = processParsed(f"../Data/ParsingResults/{p}/novor_results.tsv")
        # save novor scores
        novor_scored_df.to_csv(f"../Data/ScoringResults/{p}/novor_scored.tsv", sep='\t', index=None)

        print("Scoring ", p, " - DeepNovo")
        # score deepnovo
        deepnovo_scored_df = processParsed(f"../Data/ParsingResults/{p}/deepnovo_results.tsv")
        # save deepnovo scores
        deepnovo_scored_df.to_csv(f"../Data/ScoringResults/{p}/deepnovo_scored.tsv", sep='\t', index=None)

        print("Scoring ", p, " - DirecTag")
        # read directag
        directag_scored_df = processParsed(f"../Data/ParsingResults/{p}/direcTag_results.tsv")
        # save directag scores
        directag_scored_df.to_csv(f"../Data/ScoringResults/{p}/direcTag_scored.tsv", sep='\t', index=None)
