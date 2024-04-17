from Pipeline.AlgorithmResultParsers.AParser import AParser
import pandas as pd


class DeepNovoParser(AParser):
    def __init__(self, file: str):
        """Initializes the DeepNovoParser object.
        Parameters:
        :param file: Path to the DeepNovo result file.
        """
        self.file = file
        self.result = None

    @staticmethod
    def __cleanPeptideModifications(peptide: str) -> str:
        """Remove all modifications from a peptide sequence.
        Parameters:
        :param peptide: Peptide sequence with modifications, separated by commas.
        :return: Peptide sequence without modifications.
        """
        sequence = list()
        for aa in peptide.split(','):
            if aa.isalpha():
                sequence.append(aa[0])
        return ''.join(sequence)

    def parse(self) -> pd.DataFrame:
        """Parse the DeepNovo result file.
        :return: DataFrame containing the scan number and predicted sequence without modifications. """
        temp_df = pd.read_csv(self.file, sep='\t', header=0)
        temp_df = temp_df.drop(columns=['exact_match', 'target_seq', 'output_score', 'accuracy_AA', 'len_AA'])
        temp_df = temp_df.query('output_seq.str.contains("inf")==False')
        temp_df['output_seq'] = temp_df['output_seq'].apply(lambda x: self.__cleanPeptideModifications(x))
        self.result = pd.DataFrame({'Scan': temp_df['scan'], 'Predicted': temp_df['output_seq']})
        return self.result


if __name__ == "__main__":
    parser = DeepNovoParser(
        "../../Data/BD7_Thermo_Pool52_HCD/DeepNovo/Result_run_1.tsv")
    deepnovo_df = parser.parse()
    print(deepnovo_df)
