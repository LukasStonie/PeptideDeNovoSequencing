
class MgfPeakCounter():
    """
    A class to count the number of peaks in a .mgf file.
    """

    def __init__(self, mgf_file):
        """
        Initialize the MgfPeakCounter object.

        :param mgf_file: The path to the .mgf file to count the peaks in.
        :type mgf_file: str
        """
        self.mgf_file = mgf_file

    def count_peaks(self):
        """
        Count the number of peaks in the .mgf file.

        :return: The number of peaks in the .mgf file.
        :rtype: int
        """
        with open(self.mgf_file, 'r') as f:
            num_peaks = 0
            for line in f:
                if line.startswith('PEPMASS'):
                    num_peaks += 1
        return num_peaks