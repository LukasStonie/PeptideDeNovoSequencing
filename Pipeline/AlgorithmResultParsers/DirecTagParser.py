import pandas as pd

from Pipeline.AlgorithmResultParsers.AParser import AParser


class DirecTagParser(AParser):
    def __init__(self, file: str, skipNLines: int):
        """Initializes the DirecTagParser object.
        Parameters:
        :param file: Path to the DirecTag result file.
        :param skipNLines: Number of lines to skip at the beginning of the file. DirecTag files have a header that needs to be skipped.
        """
        self.aaList = list("ACDEFGHIKLMNPQRSTVWY")
        self.file = file
        self.skipNLines = skipNLines

    def __validSequence(self, sequence: str) -> bool:
        """Check if a sequence is valid.
        Parameters:
        :param sequence: Sequence to be checked.
        :return: True if the sequence is valid, False otherwise.
        """
        return all([aa in self.aaList for aa in sequence])
    def parse(self, remove_invalid_tags=True) -> pd.DataFrame:
        with open(self.file, 'r') as file:
            lines = file.readlines()
        entries = list()
        latestID = -1
        for idx, line in enumerate(lines):
            if idx < self.skipNLines:
                continue
            if line.startswith("S"):
                latestID = int(line.split('\t')[3])
            if line.startswith("T"):
                tagline = line.split('\t')
                entry = {'ID': latestID, 'Predicted': tagline[1], 'TotalScore': tagline[7],
                         'ComplementScore': tagline[8], 'IntensityScore': tagline[9], 'mzFidelityScore': tagline[10]}
                entries.append(entry)
        temp_df = pd.DataFrame(entries)
        if remove_invalid_tags:
            filtered_df = temp_df.apply(lambda x: self.__validSequence(x['Predicted']), axis=1)
            return temp_df[filtered_df]
        return temp_df


if __name__ == "__main__":
    parser = DirecTagParser(
        "../../Data/BD7_Thermo_Pool52_HCD/DirecTag/Run_1/01640c_BD7-Thermo_SRM_Pool_52_01_01-3xHCD-1h-R2.tags",
        25)
    novor_df = parser.parse()
    print(novor_df)
