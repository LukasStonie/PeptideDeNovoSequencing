import pandas as pd
from joblib import Parallel, delayed

import warnings

from Pipeline.Scoring.AlignmentScore import AlignmentScore
from Pipeline.Scoring.LevenshteinDistance import LevenshteinDistance
from Pipeline.Scoring.NormalizedAlignmentScore import NormalizedAlignmentScore
from Pipeline.Scoring.SequenceIdentity import SequenceIdentity
from Pipeline.Scoring.SequenceSimilarity import SequenceSimilarity

warnings.filterwarnings("ignore")

RESULT_CLOUMNS = ['ID', 'Predicted', 'Inclusion', 'Score', 'Similarity', 'Identity', 'Local Alignment', 'Global Alignment',
                  'Normalized Local Alignment', 'Normalized Global Alignment', 'Levenshtein']


def best_match(s, candidates, alignment_mode='global'):
    ''' Return the item in candidates that best matches s.

    Will return None if a good enough match is not found.
    '''
    if alignment_mode == 'local':
        similarity = SequenceSimilarity(alignment_mode=alignment_mode, open_gap_score=-4, extend_gap_score=-4)
    else:
        similarity = SequenceSimilarity(alignment_mode=alignment_mode, open_gap_score=-2, extend_gap_score=-2)
    for candidate in candidates:
        if similarity.getScore(predicted=s, actual=candidate) == 1.0:
            return [s, candidate]
    return [s, None]


def best_match_parallel(s_df, candidates, alignment_mode='global'):
    output = s_df.apply(lambda x: best_match(x['Predicted'], candidates, alignment_mode), axis=1)
    output = pd.DataFrame(list(output), columns=['Predicted', 'Inclusion'])
    return output


def calculateScores(entry, alignment_mode: str = 'global') -> list:
    identiy = SequenceIdentity(alignment_mode=alignment_mode)
    similarity = SequenceSimilarity(alignment_mode=alignment_mode)
    localAlignmentScore = AlignmentScore(alignment_mode='local')
    globalAlignmentScore = AlignmentScore(alignment_mode='global')
    normalizedLocalAlignmentScore = NormalizedAlignmentScore(alignment_mode='local')
    normalizedGlobalAlignmentScore = NormalizedAlignmentScore(alignment_mode='global')
    levenshtein = LevenshteinDistance()
    return [
        entry['ID'] if 'ID' in entry else entry['Scan'],
        entry['Predicted'], entry['Inclusion'], entry['Score'],
        similarity.getScore(predicted=entry['Predicted'], actual=entry['Inclusion']),
        identiy.getScore(predicted=entry['Predicted'], actual=entry['Inclusion']),
        localAlignmentScore.getScore(predicted=entry['Predicted'], actual=entry['Inclusion']),
        globalAlignmentScore.getScore(predicted=entry['Predicted'], actual=entry['Inclusion']),
        normalizedLocalAlignmentScore.getScore(predicted=entry['Predicted'], actual=entry['Inclusion']),
        normalizedGlobalAlignmentScore.getScore(predicted=entry['Predicted'], actual=entry['Inclusion']),
        levenshtein.getScore(predicted=entry['Predicted'], actual=entry['Inclusion'])]


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


def calculateScoresOfChunk(subset_df, alignment_mode: str = 'global'):
    output = subset_df.apply(lambda x: calculateScores(x, alignment_mode), axis=1)
    output = pd.DataFrame(list(output), columns=RESULT_CLOUMNS)
    return output


def findPerfectSimilarityMatchInclusionList(file: str, pool: str, algorithm: str, alignment_mode: str = 'global'):
    # read parsed result
    parsed_df = pd.read_csv(file, sep='\t', index_col=None, header=0).query('Actual == \' \'')
    # get unique predicted sequences and create a DataFrame
    unique_predicted_df = pd.DataFrame(parsed_df['Predicted'].unique())
    unique_predicted_df.columns = ['Predicted']

    # read the inclusion list
    inclusion_list = pd.read_csv(f"../Data/Datasets/{pool}/Thermo_SRM_{pool}_01_01_3xHCD-1h-R2-tryptic/peptides.txt",
                                 sep='\t', index_col=None, header=0)
    # reduce dfs to necessary colmns
    inclusion_list = inclusion_list[['Sequence']]

    # define number of cpus and calculate chunk indices for parallel processing
    cpus = 6  # os.cpu_count()
    chunk_indices_similarity_overlap = calculateChunkIndex(len(unique_predicted_df), cpus)

    # define parallel processing object with number of cpus
    parallel_similarity_overlap = Parallel(n_jobs=cpus)

    # process the unique predictions in parallel
    output = parallel_similarity_overlap(
        delayed(best_match_parallel)(unique_predicted_df.iloc[chunk[0]:chunk[1], :],
                                     inclusion_list['Sequence'].to_list(), alignment_mode) for chunk in
        chunk_indices_similarity_overlap)
    # concatenate the results
    output = pd.concat(output, axis=0)
    output = output.dropna()
    output.to_csv(f"../Data/ScoringResults_Unidentified/CheckInclusionList/{pool}_{algorithm}_similarity_100_match.tsv",
                  sep='\t', index=None)


def processParsed(file: str, pool: str, algorithm: str, alignment_mode: str = 'global'):
    # read parsed result
    parsed_df = pd.read_csv(file, sep='\t', index_col=None, header=0).query('Actual == \' \'')
    print(parsed_df.shape[0])
    # read the inclusion list
    sim_100_matches = pd.read_csv(
        f"../Data/ScoringResults_Unidentified/CheckInclusionList/Overlap_Similarity_100/{pool}_{algorithm}_similarity_100_match.tsv",
        sep='\t', index_col=None, header=0)
    # reduce dfs to necessary colmns

    merged_df = parsed_df.merge(sim_100_matches, left_on='Predicted', right_on='Predicted', how='left').dropna().drop(columns=['Actual'])
    print(merged_df.shape[0])
    print(merged_df.columns)

    # define number of cpus and calculate chunk indices for parallel processing
    cpus = 6  # os.cpu_count()
    chunk_indices = calculateChunkIndex(len(merged_df), cpus)

    # define parallel processing object with number of cpus
    parallel_direcTag = Parallel(n_jobs=cpus)
    # process the direcTag result in parallel
    output = parallel_direcTag(delayed(calculateScoresOfChunk)(merged_df.iloc[chunk[0]:chunk[1], :], alignment_mode) for chunk in chunk_indices)
    # concatenate the results
    output = pd.concat(output, axis=0)
    return output


if __name__ == "__main__":
    pools = ['Pool_49', 'Pool_52', 'Pool_60']
    # pools = ['Pool_49']

    if False:
        for p in pools:
            print("Finding Matches: ", p, " - PEAKS")
            findPerfectSimilarityMatchInclusionList(f"../Data/ParsingResults/{p}/peaks_results_all_sequences.tsv", p,
                                                    'peaks', alignment_mode='global')

            print("Finding Matches: ", p, " - Novor")
            findPerfectSimilarityMatchInclusionList(f"../Data/ParsingResults/{p}/novor_results_all_sequences.tsv", p,
                                                    'novor', alignment_mode='global')

            print("Finding Matches: ", p, " - DirecTag")
            findPerfectSimilarityMatchInclusionList(f"../Data/ParsingResults/{p}/direcTag_results_all_sequences.tsv", p,
                                                    'direcTag', alignment_mode='local')

            print("Finding Matches: ", p, " - DeepNovo")
            findPerfectSimilarityMatchInclusionList(f"../Data/ParsingResults/{p}/deepnovo_results_all_sequences.tsv", p,
                                                    'deepnovo', alignment_mode='local')

    for p in pools:
        print("Scoring ", p, " - PEAKS")
        # score peaks
        peaks_scored_df = processParsed(f"../Data/ParsingResults/{p}/peaks_results_all_sequences.tsv", p, 'peaks',alignment_mode='global')
        # save peaks scores
        peaks_scored_df.to_csv(f"../Data/ScoringResults_Unidentified/CheckInclusionList/{p}/peaks_scored.tsv", sep='\t', index=None)

        print("Scoring ", p, " - Novor")
        # score novor
        novor_scored_df = processParsed(f"../Data/ParsingResults/{p}/novor_results_all_sequences.tsv",p, 'novor', alignment_mode='global')
        # save novor scores
        novor_scored_df.to_csv(f"../Data/ScoringResults_Unidentified/CheckInclusionList/{p}/novor_scored.tsv", sep='\t',index=None)

        print("Scoring ", p, " - DeepNovo")
        # score deepnovo
        #deepnovo_scored_df = processParsed(f"../Data/ParsingResults/{p}/deepnovo_results_all_sequences.tsv",p,'deepnovo', alignment_mode='local')
        # save deepnovo scores
        #deepnovo_scored_df.to_csv(f"../Data/ScoringResults_Unidentified/CheckInclusionList/{p}/deepnovo_scored.tsv", sep='\t', index=None)

        print("Scoring ", p, " - DirecTag")
        # read directag
        #directag_scored_df = processParsed(f"../Data/ParsingResults/{p}/direcTag_results_all_sequences.tsv",p,'direcTag', alignment_mode='local')
        # save directag scores
        #directag_scored_df.to_csv(f"../Data/ScoringResults_Unidentified/CheckInclusionList/{p}/direcTag_scored.tsv", sep='\t', index=None)
