

class SAM(object):

    def __init__(
            self,
            file_name):
        self.file_name = file_name
        self.file_handle = open(self.file_name)
        self.header = []
        self.targets = []

    def _parse_header(self):
        for line in self.file_handle:
            if line[0] == "@":
                # Its a header line
                fields = line.split()
                if fields[0] == "@SQ":
                    self.targets[fields[1]] = len(self.targets), fields[2:]
                else:
                    self.header.append(fields[1:])
            else:
                # Reset file handle to start of first alignment record
                pos = self.file_handle.tell()
                pos -= len(line)
                self.file_handle.seek(pos)

    def next(self):
        # Iterate through alignment



