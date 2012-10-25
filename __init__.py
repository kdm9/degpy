from collections import deque
import sqlite3

#################################Hard Config##################################
CIGAR_MAP = {
    0: "M",
    1: "I",
    2: "D",
    3: "N",
    4: "S",
    5: "H",
    6: "P",
    7: "=",
    8: "X"
    }

# Minimum number of matches at the 5' end of a read's mapping position
MIN_INITIAL_MATCH = 5

# Minimum MAPQ value, see SAM documentation:
# MAPQ =  -10 log10 Pr{mapping position is wrong}, rounded to the nearest
# integer. A value 255 indicates that the mapping quality is not available.
MIN_MAP_QUALITY = -1

# QC data: this will be checked against the tags in SAM file below
MAX_MISMATCHES = 1  # XM
MAX_GAP_OPENS = 0  # XO
MAX_GAP_EXTENSIONS = 0  # XG
MIN_SCORE = 20  # AS

#how many points across gene length
BINS = 100.0

# Minimum number of hits per transcript before a histogram is made
MIN_HITS_PER_TX = int(BINS / 2.0)

# Roughly MAX of 500Mb assuming entry size of 100b
MAX_LINES_IN_MEM = 5000000


def tuple_list_to_dict(this_list):
    this_dict = {}
    for item in this_list:
        key, value = item[0], item[1:]
        this_dict[key] = value
    return this_dict


def dict_to_tuple_list(this_dict):
    this_tuple_list = []
    for key, value in this_dict:
        if not isinstance(value, tuple):
            value = (value,)
        this_tuple_list.append((key,) + value)
    return this_tuple_list


def average(data):
    data = list(data)
    return float(sum(data)) / float(len(data))


def moving_average(data, window_size=5):
    data = list(data)
    for i in xrange(window_size):
        data.append(0.0)
    moving_avg = []
    window = deque([0.0] * window_size)
    o = int(window_size / 2)
    for value in data:
        window.append(value)
        window.popleft()
        moving_avg.append(average(window))
    return moving_avg[o:-(o + 1)]


def stdev(data):
    data = list(data)
    avg = average(data)
    sum_squares = 0.0
    for x in data:
        sum_squares += ((float(x) - avg) ** 2)
    return sum_squares / float(len(data) - 1)


def moving_stdev(data, window_size=5):
    data = list(data)
    for i in xrange(window_size):
        data.append(0.0)
    moving_sd = []
    window = deque([0.0] * window_size)
    o = int(window_size / 2)  # Offset, remove the first o and last o+1 values
    for value in data:
        window.append(value)
        window.popleft()
        moving_sd.append(stdev(window))
    return moving_sd[o:-(o + (window_size - 2 * o))]  # allows even window_size


def qc_hit(read_mapq, cigar_list, tags):
    # If the map does not start with a match
    tag_dict = tuple_list_to_dict(tags)
    if cigar_list[0][0] != 0:
        return False
    # If there are not at least MIN_INITIAL_MATCH matches at the 5' end of
    # of the map
    elif cigar_list[0][1] < MIN_INITIAL_MATCH:
        return False
    elif tag_dict['XO'] < MAX_GAP_OPENS:
        return False
    elif tag_dict['XG'] < MAX_GAP_EXTENSIONS:
        return False
    elif tag_dict['XM'] < MAX_MISMATCHES:
        return False
    #elif tag_dict['MD']:
    #    pass
    else:
        return True


def setup_database():
    pass


def make_cigar_string(cigar_list):
        cigar_str = ""
        for section in cigar_list:
            # section is tuple (cigar_type, count)
            cigar_str += str(section[1]) + CIGAR_MAP[section[0]]
        return cigar_str


class Base(object):

    def _init_db(self):
        self.sqlite_con = sqlite3.connect(self.db_file_name)

    def _setup_db(self):
        ## Write tables to file db

        cur = self.sqlite_con.cursor()

        # samples table, everything is tied to an sample
        create_samples_table = """CREATE TABLE IF NOT EXISTS samples
        (id INTEGER PRIMARY KEY AUTOINCREMENT, sam_file_name TEXT, date TEXT,
        sample_type INTEGER NOT NULL)"""
        cur.execute(create_samples_table)
        self.sqlite_con.commit()

        # Targets table, holds mapping between AGI and tid from SAM file header
        create_target_table = """CREATE TABLE IF NOT EXISTS targets
        (id INTEGER, agi TEXT, length INTEGER, hits INTEGER,
        sample_id INTEGER, bin_width REAL, FOREIGN KEY(sample_id)
        REFERENCES analyses(id))"""
        cur.execute(create_target_table)
        self.sqlite_con.commit()

        # Hits table, holds mappings between reads and AGIs
        create_hits_table = """CREATE TABLE IF NOT EXISTS hits (read_name TEXT,
        agi_id INTEGER, pos INTEGER, cigar TEXT, mapq INTEGER,
        passed_qc INTEGER, sample_id INTEGER,  FOREIGN KEY(sample_id)
        REFERENCES analyses(id), FOREIGN KEY(agi_id) REFERENCES targets(id))"""
        cur.execute(create_hits_table)
        self.sqlite_con.commit()

        # Histograms table
        create_histograms_table = """CREATE TABLE IF NOT EXISTS histograms
        (id INTEGER PRIMARY KEY, hist_str TEXT, peak_str TEXT, length INTEGER,
        sample_id INTEGER, FOREIGN KEY (sample_id) REFERENCES
        analyses(id), FOREIGN KEY(agi_id) REFERENCES targets(id))"""
        cur.execute(create_histograms_table)
        self.sqlite_con.commit()

        # Define SQL statements
        ## samples
        self.insert_sample = """INSERT INTO samples (sam_file_name, date,
        sample_type) VALUES (?, datetime('now'), ?)"""
        self.select_sample_id = """SELECT id FROM samples WHERE sam_file_name
        = ?"""
        self.select_sam_file = """SELECT sam_file_name FROM samples WHERE id
        = ?"""

        ## targets
        self.insert_new_target = """INSERT INTO targets (id, agi, length, hits,
        sample_id, bin_width) VALUES (?, ?, ?, 0, ?, ?)"""
        self.select_all_tids = """SELECT id FROM targets WHERE sample_id = ?"""
        self.select_targets = """SELECT agi, length FROM targets WHERE id
        = ? AND sample_id = ?"""
        self.select_target_lengths = """SELECT id, length FROM targets WHERE
        sample_id = ?"""
        self.update_target_hits = """UPDATE targets SET hits = ? WHERE id =
         ? AND sample_id = ?"""
        self.select_targets_by_hits = """SELECT id FROM targets WHERE sample_id
        = ? ORDER BY hits DESC LIMIT ?"""

        ## hits
        self.insert_hit = """INSERT INTO hits (read_name, agi_id, pos, cigar,
        mapq, passed_qc, sample_id) VALUES (?, ?, ?, ?, ?, ?, ?)"""
        self.select_hits_for_target = """SELECT read_name, agi_id, pos, cigar,
        mapq, passed_qc FROM hits ORDER BY pos ASC"""
        self.select_all_hits = """SELECT agi_id, pos, read_name, cigar, mapq,
        passed_qc FROM hits WHERE sample_id = ? ORDER BY agi_id ASC, pos ASC"""
        self.select_hit_read_names = """SELECT read_name FROM hits WHERE
        sample_id = ?"""

        ## histograms
        self.insert_histogram = """INSERT INTO histograms (agi_id, hist_str,
        length, peak_str, sample_id) VALUES (?, ?, ?, ?, ?)"""
        self.select_all_histograms = """SELECT hist_str, length, peak_str FROM
        histograms WHERE sample_id = ?"""
        self.select_agi_bins = """SELECT bin_num, bin_count FROM histograms
        WHERE sample_id = ? AND agi_id = ?"""

        cur.close()
