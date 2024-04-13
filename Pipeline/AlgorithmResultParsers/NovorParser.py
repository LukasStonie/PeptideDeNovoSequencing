from Pipeline.Scoring import SequenceIdentity, SequenceSimilarity
import pandas as pd

class NovorParser():
    def __init__(self, file:str, skipNLines:int):
        """Initializes the NovorParser object.
        Parameters:
        :param file: Path to the Novor result file.
        :param skipNLines: Number of lines to skip at the beginning of the file. Novor files have a header that needs to be skipped.
        :param scanNumberCorrection: Number to be added to the scan number. Novor uses scan number 0 for all scans, so this number is used to correct it.

        """
        self.file = file
        self.skipNLines = skipNLines

        self.novor_df = None

        self.__parse()

    def __parse(self):
        temp_df = pd.read_csv(self.file, skiprows=self.skipNLines, sep=',', header=0)
        temp_df.columns = temp_df.columns.str.strip('# ')
        temp_df['id']= temp_df['id']-1
        temp_df=temp_df.drop(columns=[''])
        self.novor_df = temp_df

    def getParsedData(self):
        return self.novor_df


def calculateScores(entry):
    identiy = SequenceIdentity()
    similarity = SequenceSimilarity()
    print(entry['peptide'], entry['Sequence'], similarity.getScore(entry['peptide'], entry['Sequence']), identiy.getScore(entry['peptide'], entry['Sequence']) )

if __name__ == "__main__":
    novor = NovorParser("/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/BD7_Thermo_Pool52_HCD/Novor/Run_2/01640c_BD7-Thermo_SRM_Pool_52_01_01-3xHCD-1h-R2.mzml.novor.csv", 20)
    expected = pd.read_csv("/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/BD7_Thermo_Pool52_HCD/Thermo_SRM_Pool_52_01_01_3xHCD-1h-R2-tryptic/msmsScans.txt", sep='\t')
    expected_seq = pd.DataFrame(expected['Sequence'])
    merged_df = pd.merge(novor.getParsedData(), expected_seq, left_on='id', right_index=True, how='inner')
    print(merged_df)
    merged_df.head(100).apply(lambda x: calculateScores(x), axis=1)
