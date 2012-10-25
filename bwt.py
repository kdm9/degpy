from array import array
# Cheap-arse FM-index

class BWT(object):
    """
    Performs the Burrows-Wheeler transform to map small reads to a reference.
    Firstly, instantiate the class. Then explictly call transform(). Then
    locate() can be called many times to find the index of a substring. read()
    and write() read and write the transform to a file to save repeated
    computation.
    """

    def __init__(self, sequence):
        self.seq = sequence
        self.seq_len = len(self.seq)
        pass

    def transform(self):
        pass

    def locate(self, substring):
        pass

    def read(self, filename):
        pass

    def write(self, filename):
        pass

from array import array

s = "ajvkdabvakgjabna"


def burrows_string(string, index):
    l = len(string)
    return string[index:l] + string[0:index]


def indxa(string):
    string += '\0'
    strlen = len(string)
    # Make-wrap-around
    S = []
    C = {}
    for iii in xrange(strlen):
        f = string[iii]
        l = string[(iii - strlen) % strlen - 1]
        S.append((f, iii, l))
    S.sort(key=lambda a: a[0])
    F = []
    last_f = None
    ctr = 0
    for f, i, l in S:
        if f != last_f:
            C[f] = ctr
        print f, burrows_string(string, i)
        F.append(f)
        last_f = f
        ctr += 1
    return S, C

indxa(s)


def finda(substr, S):
    rss = substr[::-1] #rev substr
    alpha_indicies = []

