from Pipeline.Scoring import SequenceIdentity, SequenceSimilarity
import pandas as pd

class DeepNovoParser():
    def __init__(self, file:str):
        """Initializes the DeepNovoParser object.
        Parameters:
        :param file: Path to the DeepNovo result file.
        """
        self.file = file

        self.deepnovo_df = None

        self.__parse()


    @staticmethod
    def __removeModifications(self, aa_list:str)->str:
        sequence = list()
        for aa in aa_list.split(','):
            if aa.isalpha():
                sequence.append(aa[0])
        return ''.join(sequence)
    def __parse(self):
        temp_df = pd.read_csv(self.file, sep='\t', header=0)
        temp_df=temp_df.drop(columns=['exact_match'])
        temp_df = temp_df.query('output_seq.str.contains("inf")==False')
        temp_df['output_seq'] = temp_df['output_seq'].apply(lambda x: ''.join([aa[0] for aa in x.split(',') if aa.isalpha()]))
        self.deepnovo_df = temp_df

    def getParsedData(self):
        return self.deepnovo_df


def calculateScores(entry):
    identiy = SequenceIdentity()
    similarity = SequenceSimilarity()
    output_seq = ''.join([aa[0] for aa in entry['output_seq'] if aa.isalpha()])
    print(output_seq, entry['Sequence'], similarity.getScore(output_seq, entry['Sequence']), identiy.getScore(output_seq, entry['Sequence']) )

if __name__ == "__main__":
    deepnovo_df = DeepNovoParser("/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/Result/Analysis_DeepNovo.tsv")
    expected = pd.read_csv("/Users/lukas/University/Bachelor_Thesis/Project/PeptideDeNovoSequencing/Data/BD7_Thermo_Pool52_HCD/Thermo_SRM_Pool_52_01_01_3xHCD-1h-R2-tryptic/msmsScans.txt", sep='\t')
    expected_seq = pd.DataFrame({'Scan number': expected['Scan number'], 'Sequence': expected['Sequence']})
    merged_df = pd.merge(deepnovo_df.getParsedData(), expected_seq, left_on='scan', right_on='Scan number', how='inner')
    print(merged_df['output_seq'])
    #res = merged_df.query('Sequence != \' \'').apply(lambda x: calculateScores(x), axis=1)
    #analysis = pd.DataFrame(list(res), columns=['Peptide', 'Sequence', 'Similarity', 'Identity'])
    #print(analysis)