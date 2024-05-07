import pandas as pd

from Pipeline.AlgorithmResultParsers.AParser import AParser


class DirecTagParser(AParser):
    def __init__(self, file: str, skipNLines: int):
        """Initializes the DirecTagParser object.
        Parameters:
        :param file: Path to the DirecTag result file.
        :param skipNLines: Number of lines to skip at the beginning of the file. DirecTag files have a header that needs to be skipped.
        """
        super().__init__(file)
        self.skipNLines = skipNLines
        self.modification_mapping = dict()

    def __validSequence(self, sequence: str) -> bool:
        """Check if a sequence is valid.
        Parameters:
        :param sequence: Sequence to be checked.
        :return: True if the sequence is valid, False otherwise.
        """
        return all([aa in self.aaList for aa in sequence])

    @staticmethod
    def __cleanPeptideModifications(peptide: str, modification_mapping:dict) -> str:
        """Remove all modifications from a peptide sequence.
        Parameters:
        :param peptide: Peptide sequence with modifications.
        :param modification_mapping: Mapping from a number to an amino acid.
        :return: Peptide sequence without modifications.
        """
        return ''.join([modification_mapping[int(aa)] if aa.isdigit() else aa for aa in peptide])
    def parseVariableModifications(self, line: str) -> dict:
        """Parse the variable modifications mapping from the result file.
        Parameters:
        :param line: Line of result file containing variable modifications mappings.
        :return: Mapping from a number to an amino acid.
        """
        # modifcations are in the format:
        # DynamicMods: D 0 37.946941 M 1 15.994915

        # split the line by ': ' and get the second part
        modifications = line.split(', DynamicMods: ')[1].split(',')[0]
        # split the modifications by ' ' and create a list of 3 elements
        modifications = modifications.split(' ')
        modifications = [modifications[i:i + 3] for i in range(0, len(modifications), 3)]
        # create a dictionary mapping the number to the amino acid
        mapping = dict()
        for mod in modifications:
            mapping[int(mod[1])] = mod[0]
        return mapping
    def parse(self,) -> pd.DataFrame:
        with open(self.file, 'r') as file:
            lines = file.readlines()
        entries = list()
        latestID = -1
        for idx, line in enumerate(lines):
            if ", DynamicMods: " in line:
                self.modification_mapping = self.parseVariableModifications(line)
                continue
            if idx < self.skipNLines:
                continue
            if line.startswith("S"):
                latestID = int(line.split('\t')[3])
            if line.startswith("T"):
                tagline = line.split('\t')
                entry = {'ID': latestID, 'Predicted': tagline[1], 'Score': tagline[7],
                         'ComplementScore': tagline[8], 'IntensityScore': tagline[9], 'mzFidelityScore': tagline[10]}
                entries.append(entry)
        temp_df = pd.DataFrame(entries)
        temp_df['Predicted'] = temp_df.apply(lambda x: self.__cleanPeptideModifications(x['Predicted'], self.modification_mapping), axis=1)
        self.result = temp_df
        return self.result


if __name__ == "__main__":
    parser = DirecTagParser(
        "../../Data/AlgorithmResults/Pool_49/DirecTag/Run_1/01640c_BA7-Thermo_SRM_Pool_49_01_01-3xHCD-1h-R2.tags",
        25)
    novor_df = parser.parse()
    print(novor_df)
