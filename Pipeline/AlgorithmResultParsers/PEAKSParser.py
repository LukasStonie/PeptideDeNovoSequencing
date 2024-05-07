import pandas as pd

from Pipeline.AlgorithmResultParsers.AParser import AParser
from Pipeline.AlgorithmResultParsers.IAverageAminoAcidScore import IAverageAminoAcidScore


class PEAKSParser(AParser, IAverageAminoAcidScore):
    def __init__(self, file: str):
        """Initializes the PEAKSParser object.
        Parameters:
        :param file: Path to the PEAKS result file.
        """
        super().__init__(file)
        self.result = None

    @staticmethod
    def getAverageAAScore(aa_list: list) -> float:
        """Calculate the average amino acid score of a list of amino acids.
        Parameters:
        :param aa_list: List of amino acids.
        :return: Average amino acid score.
        """
        return sum(aa_list) / len(aa_list)

    @staticmethod
    def __cleanPeptideModifications(peptide: str) -> str:
        """Remove all modifications from a peptide sequence.
        Parameters:
        :param peptide: Peptide sequence with modifications.
        :return: Peptide sequence without modifications.
        """
        return ''.join([aa for aa in peptide if aa.isalpha()])

    def parse(self) -> pd.DataFrame:
        temp_df = pd.read_csv(self.file, header=0, index_col=False)
        temp_df['Peptide'] = temp_df['Peptide'].apply(lambda x: self.__cleanPeptideModifications(x))
        temp_df['local confidence (%)'] = temp_df['local confidence (%)'].apply(lambda x: self.getAverageAAScore([int(score) for score in x.split(' ')]))
        self.result = pd.DataFrame({'Scan': temp_df['Scan'], 'Predicted': temp_df['Peptide'], 'Score': temp_df['local confidence (%)']})
        return self.result


if __name__ == "__main__":
    parser = PEAKSParser("../../Data/BD7_Thermo_Pool52_HCD/PEAKS/AlgorithmResults 1/Sample 1.denovo.csv")
    peaks_df = parser.parse()
    print(peaks_df)
