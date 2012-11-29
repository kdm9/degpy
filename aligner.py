from subprocess import Popen, PIPE
import subprocess
from degpy import exe_path, ExeNotInstalled


class Aligner(object):
    
    def __init__(self, exes):
        self.exepaths = []
        for exe in exes:
            exepath = exe_path(exe)
            if exepath is None:
                raise RuntimeError("Not installed: %s" % exe)
            self.exepaths.append(exepath)


class BWA(Aligner):
    """
    """
    exes = ["bwa", "samtools"]  # 
    def __init__(self):
        super(BWA, self).__init__(exes)

    def run(self):
        pass

