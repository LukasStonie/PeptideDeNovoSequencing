from Pipeline.AlgorithmResultParsers.AParser import AParser
import pandas as pd


class NovorParser(AParser):
    def __init__(self, file: str, skipNLines: int):
        """Initializes the NovorParser object.
        Parameters:
        :param file: Path to the Novor result file.
        :param skipNLines: Number of lines to skip at the beginning of the file. Novor files have a header that needs to be skipped.
        :param scanNumberCorrection: Number to be added to the scan number. Novor uses scan number 0 for all scans, so this number is used to correct it.

        """
        self.file = file
        self.skipNLines = skipNLines
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
        """Parse the Novor result file.
        :return: DataFrame containing the scan number, predicted sequence without modifications.
        """
        temp_df = pd.read_csv(self.file, skiprows=self.skipNLines, sep=',', header=0)
        temp_df.columns = temp_df.columns.str.strip('# ')
        temp_df['id'] = temp_df['id'] - 1
        temp_df = temp_df.drop(
            columns=['', 'RT', 'score', 'aaScore', 'ppm(1e6*err/(mz*z))', 'err(data-denovo)', 'pepMass(denovo)', 'z',
                     'mz(data)', 'scanNum'])
        temp_df['peptide'] = temp_df['peptide'].apply(lambda x: self.__cleanPeptideModifications(x))
        self.result = pd.DataFrame({'ID': temp_df['id'], 'Predicted': temp_df['peptide']})
        return self.result


if __name__ == "__main__":
    parser = NovorParser(
        "../../Data/BD7_Thermo_Pool52_HCD/Novor/Run_2/01640c_BD7-Thermo_SRM_Pool_52_01_01-3xHCD-1h-R2.mzml.novor.csv",
        20)
    novor_df = parser.parse()
    print(novor_df)
