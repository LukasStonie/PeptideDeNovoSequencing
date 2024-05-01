import pandas as pd

from Pipeline.AlgorithmResultParsers.AParser import AParser


class PEAKSParser(AParser):
    def __init__(self, file: str):
        """Initializes the PEAKSParser object.
        Parameters:
        :param file: Path to the PEAKS result file.
        """
        self.file = file
        self.result = None

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
        self.result = pd.DataFrame({'Scan': temp_df['Scan'], 'Predicted': temp_df['Peptide']})
        return self.result


if __name__ == "__main__":
    parser = PEAKSParser("../../Data/BD7_Thermo_Pool52_HCD/PEAKS/AlgorithmResults 1/Sample 1.denovo.csv")
    peaks_df = parser.parse()
    print(peaks_df)
