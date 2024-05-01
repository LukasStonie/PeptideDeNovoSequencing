from typing import Tuple, List


class MGFEntry:
    def __init__(self, sequence:str):
        self.title = 'title'
        self.sequence = sequence
        self.scan = -1
        self.charge = -1
        self.rtinseconds = -1
        self.pepmass = -1
        self.peaks = list()

    def setTitle(self, title:str):
        self.title = title

    def setSequence(self, sequence:str):
        self.sequence = sequence

    def setScan(self, scan:int):
        self.scan = scan

    def setCharge(self, charge:int):
        self.charge = charge

    def setRTInSeconds(self, rtinseconds:float):
        self.rtinseconds = rtinseconds

    def setPepMass(self, pepmass:float):
        self.pepmass = pepmass

    def addPeak(self, peak:Tuple[float, float]):
        self.peaks.append(peak)

    def __str__(self):
        peaks = '\n'.join([f"{peak[0]} {peak[1]}" for peak in self.peaks])
        return f"BEGIN IONS\nTITLE={self.title}\nPEPMASS={self.pepmass}\nCHARGE={self.charge}+\nSCANS={self.scan}\nRTINSECONDS={self.rtinseconds}\nSEQ={self.sequence}\n{peaks}\nEND IONS\n"

class DeepNovoPreProcessor:
    def __init__(self, destination:str, sequence:str):
        self.destination = destination
        self.sequence = sequence

    def process(self, source:str):
        entries = self.__parse(source)
        self.__write(entries)

    def __parse(self, source:str)->List[MGFEntry]:
        lines = list()
        with open(source, 'r') as file:
            lines = file.readlines()
        entries = list()
        entry = None
        for line in lines:
            if line.startswith("BEGIN IONS"):
                entry = MGFEntry(self.sequence)
            elif line.startswith("TITLE"):
                entry.setTitle(line.split('=',1)[1].strip())
            elif line.startswith("PEPMASS"):
                entry.setPepMass(float(line.split('=')[1].strip()))
            elif line.startswith("CHARGE"):
                entry.setCharge(int(line.split('=')[1].strip().split('+')[0]))
            elif line.startswith("SCANS"):
                entry.setScan(int(line.split('=')[1].strip()))
            elif line.startswith("RTINSECONDS"):
                entry.setRTInSeconds(float(line.split('=')[1].strip()))
            elif line.startswith("END IONS"):
                entries.append(entry)
            else:
                mz, intensity =line.split()
                entry.addPeak(tuple([float(mz), float(intensity)]) )
        return entries

    def __write(self, entries):
        with open(self.destination, 'w') as file:
            for entry in entries:
                file.write(str(entry))


if __name__ == "__main__":
    files = ['Pool_49/01640c_BA7-Thermo_SRM_Pool_49_01_01-3xHCD-1h-R2', 'Pool_52/01640c_BD7-Thermo_SRM_Pool_52_01_01-3xHCD-1h-R2','Pool_60/01640c_BD8-Thermo_SRM_Pool_60_01_01-3xHCD-1h-R2']
    for file in files:
        file = f"../../Data/Datasets/{file}"
        preprocessor = DeepNovoPreProcessor(file+'_deepnovo.mgf', "PEPTIDE")
        preprocessor.process(file+'.mgf')