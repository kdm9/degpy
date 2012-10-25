from matplotlib import pyplot
import DegradomeAnalyser as da


def make_histogram(cursor):
    histogram = [0 for iii in xrange(int(da.BINS))]
    these_bins = cursor.fetchmany()
    while these_bins:
        for bin_tuple in these_bins:
            (bin_num, bin_count) = bin_tuple
            histogram[bin_num] += bin_count
        these_bins = cursor.fetchmany()
    return histogram


class Analysis(object):

    def __init__(self, db_file_name, output_dir, verbose=False):
        self.db_file_name = db_file_name
        self.output_dir = output_dir
        self.verbose = verbose
        self.counts_dict = {}
        self._init_db()

    def _summarise_sample_degradome(self, sample_id):
        if self.verbose:
            print ""
        cur = self.sqlite_con.cursor()
        cur.arraysize = da.MAX_LINES_IN_MEM
#        peaks = [0 for iii in xrange(int(da.BINS))]
        histogram_cur = cur.execute(
            self.select_all_bins,
            (sample_id,)
            )
        histogram = make_histogram(histogram_cur)
        self.counts_dict[sample_id] = sum(histogram)

        cur.close()

    def make_culmulative_histogram(self, sample_id):
        if self.verbose:
            print "make_culmulative_histogram"
        cur = self.sqlite_con.cursor()
        cur.arraysize = da.MAX_LINES_IN_MEM
        histogram = [0 for iii in xrange(int(da.BINS))]
        culm_histogram = []
        histogram_cur = cur.execute(
            self.select_all_bins,
            (sample_id,)
            )
        these_bins = histogram_cur.fetchmany()
        while these_bins:
            for bin_tuple in these_bins:
                (bin_num, bin_count) = bin_tuple
                histogram[bin_num] += bin_count
            these_bins = histogram_cur.fetchmany()
        culm_sum = 0
        for iii in xrange(len(histogram)):
            culm_sum += histogram[iii]
            culm_histogram.append(culm_sum)
        print culm_histogram
        pyplot.clf()
        pyplot.plot(range(len(culm_histogram)), culm_histogram)
        pyplot.suptitle(str("Overall culmulative"))
        pyplot.savefig(self.output_dir + "%i_overall_culm.png" % sample_id)
        cur.close()

    def make_individual_histograms(self, sample_id):
        cur = self.sqlite_con.cursor()
        cur.arraysize = da.MAX_LINES_IN_MEM
        targets_cur = cur.execute(
            self.select_targets_by_hits,
            (sample_id, 100)
            )
        these_targets = targets_cur.fetchmany()
        while these_targets:
            for (tid,) in these_targets:
                self._get_histogram_for_agi_id(tid, sample_id)
            these_targets = targets_cur.fetchmany()

    def make_tplot(self, sample_id, target_id):
        cur = self.sqlite_con.cursor()
        cur.arraysize = da.MAX_LINES_IN_MEM
#        histogram = [0.0 for iii in xrange(100)]

    def normalise_histogram(self, histogram):
        pass

    def _pick_peaks(self, histogram):
        length = len(histogram)
        window = 5
        for iii in xrange(length - window):
            pass

    def pick_all_peaks(self):
        pass

    def _get_histogram_for_agi_id(self, agi_id, sample_id):
        cur = self.sqlite_con.cursor()
        histogram_cur = cur.execute(
            self.select_agi_bins,
            (sample_id, agi_id)
            )
        histogram = [0 for iii in xrange(int(da.BINS))]
        these_bins = histogram_cur.fetchmany(da.MAX_LINES_IN_MEM)
        while these_bins:
            for bin_tuple in these_bins:
                (bin_num, bin_count) = bin_tuple
                histogram[bin_num] += bin_count
            these_bins = histogram_cur.fetchmany(da.MAX_LINES_IN_MEM)
        print agi_id, "\n", histogram
        pyplot.clf()
        pyplot.plot(range(len(histogram)), histogram)
        pyplot.suptitle(str(agi_id))
        pyplot.savefig(self.output_dir + "%i_%s.png" % (sample_id, agi_id))
        cur.close()
