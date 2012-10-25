import pysam
import DegradomeAnalyser as da


class Loader(object):
    """
    Degradome profile examiner
    To run, first initialise
    """
    def __init__(self, db_file_name, output_dir, verbose=False):
        self.db_file_name = db_file_name
        self.output_dir = output_dir
        self.verbose = verbose
        self._setup_db()

    def _get_sam_file_name(self, sample_id):
        cur = self.sqlite_con.cursor()
        sam_file_name = cur.execute(
            self.select_sam_file,
            (sample_id,)
            ).fetchone()[0]
        return sam_file_name

    def _get_sam_file_obj(self, sample_id, mode="r"):
        sam_file_name = self._get_sam_file_name(sample_id)
        return pysam.Samfile(sam_file_name, mode)

    def add_sample(self, sam_file_name, sample_type):
        """
        Checks if sam file is in database, if it is not, it is added, and
        in either case the sample id is returned
        """
        cur = self.sqlite_con.cursor()
        sample_id = cur.execute(
            self.select_sample_id,
            (sam_file_name,)
            ).fetchone()
        if sample_id is None:
            cur.execute(self.insert_sample, (sam_file_name, sample_type))
            self.sqlite_con.commit()
            sample_id = cur.execute("SELECT last_insert_rowid()").fetchone()[0]
            self.sqlite_con.commit()
        else:
            sample_id = sample_id[0]
        cur.close()
        return sample_id

    def fill_targets_table(self, sample_id):
        cur = self.sqlite_con.cursor()
        cur.arraysize = da.MAX_LINES_IN_MEM
        sam_file = self._get_sam_file_obj(sample_id)
        target_agis = sam_file.references
        target_lengths = sam_file.lengths

        # Get exsiting tids
        tid_cur = cur.execute(
            self.select_all_tids,
            (sample_id,)
            ).fetchall()
        existing_tids = []
        for tid in tid_cur:
            existing_tids.append(tid[0])

        targets = []
        for ttt in xrange(len(target_agis)):
            target_agi = target_agis[ttt]
            tid = sam_file.gettid(target_agi)
            if tid not in existing_tids:
                target_length = target_lengths[ttt]
                bin_width = float(target_length) / da.BINS
                targets.append(
                    (tid, target_agi, target_length, sample_id, bin_width)
                    )
        if len(targets) > 0:
            cur.executemany(self.insert_new_target, targets)
            self.sqlite_con.commit()
        if self.verbose:
            print "wrote %i targets" % len(targets)
            print "last was: ", targets[-1:]
        cur.close()
        sam_file.close()

    def fill_hits_table(self, sample_id):
        """
        """
        cur = self.sqlite_con.cursor()
        cur.arraysize = da.MAX_LINES_IN_MEM
        hits = []
        sam_file = self._get_sam_file_obj(sample_id)
        # Get existing reads. There shouldn't be any, but its good to check
        read_name_cour = cur.execute(self.select_hit_read_names, (sample_id,))
        existing_read_names = read_name_cur.fetchall()
        read_count = 0
        for read in sam_file:
            read_name = read.qname
            if read_name in existing_read_names:
                continue
            if not read.is_unmapped:
                read_mapq = read.mapq
                read_cigar = read.cigar
                read_cigar_str = self._make_cigar_string(read_cigar)
                read_passed_qc = self._qc_map(read_mapq, read_cigar, read.tags)
                hit = (
                       read_name,
                       read.tid,
                       read.pos,
                       read_cigar_str,
                       read_mapq,
                       read_passed_qc,
                       sample_id
                      )
                hits.append(hit)
            # Manage Memory, keeps mem usage below MAX_LINES_IN_MEM
            if len(hits) > da.MAX_LINES_IN_MEM:
                read_count += len(hits)
                if self.verbose:
                    print "inserting %i mapped reads, buffer is full" % \
                     (len(hits), )
                cur.executemany(self.insert_hit, hits)
                self.sqlite_con.commit()
                del hits[:]
        read_count += len(hits)
        cur.executemany(self.insert_hit, hits)
        self.sqlite_con.commit()
        cur.close()
        if self.verbose:
            print "inserted %i mapped reads" % (len(hits), )
            print "got %i hits" % (read_count, )
        sam_file.close()

    def remove_histogram_peaks(self, histogram):
        pass

    def fill_histograms_table(self, sample_id):
        """
        """
        updated_targets = []
        histograms = []
        bin_width = (float(da.BINS) / 100.0)
        print "bin withd is %f " % bin_width

        # Get target lengths and convert to dict
        cur = self.sqlite_con.cursor()
        cur.arraysize = da.MAX_LINES_IN_MEM
        targets = cur.execute(
                              self.select_target_lengths,
                              (sample_id,)
                             ).fetchall()
        targets = self._tuple_list_to_dict(targets)

        # Get hits from DB
        hits_cur = cur.execute(
                               self.select_all_hits,
                                   (sample_id,)
                              )

        # one chunk at a time
        these_hits = hits_cur.fetchmany()
        # Initialise per-agi variables.
        last_agi_id = -1
        norm_hit_plot = []
        binned_hit_plot = [0]
        bin_num = 0
        while these_hits:
            for hhh in xrange(len(these_hits)):
                hit = these_hits[hhh]
                (agi_id, pos, read_name, cigar, mapq, passed_qc) = hit
                # Do this once per agi
                if agi_id > last_agi_id:
                    # If this is not the first read
                    if last_agi_id > -1:
                        # exclude targets with too few hits
                        if len(norm_hit_plot) > da.MIN_HITS_PER_TX:
                            # Update targets table with correct data, format is
                            # (hits, agi_id, sample_id)
                            this_agi_hits = len(norm_hit_plot)
                            updated_targets.append(
                                (this_agi_hits, agi_id, sample_id)
                                )
                            # Fill histogram
                            next_barrier = bin_width
                            for norm_hit in norm_hit_plot:
                                # find the bin of the current read
                                while norm_hit > next_barrier \
                                 and bin_num < da.BINS:
                                        bin_num += 1
                                        next_barrier += bin_width
                                        binned_hit_plot.append(0)
                                # and increment it
                                binned_hit_plot[bin_num] += 1
                            # Fill histogram to full length
                            while len(binned_hit_plot) < int(da.BINS):
                                binned_hit_plot.append(0)
                            # Append histogram to table
                            for b in xrange(len(binned_hit_plot)):
                                histograms.append(
                                    (agi_id, b, binned_hit_plot[b], sample_id)
                                    )
                            # Keep memory usage low.
                            if len(histograms) > da.MAX_LINES_IN_MEM:
                                cur.executemany(
                                    self.insert_histogram,
                                    histograms
                                    )
                                self.sqlite_con.commit()
                                del histograms[:]
                    # setup new values
                    (length,) = targets[agi_id]
                    bin_num = 0
                    binned_hit_plot = [0]
                    del norm_hit_plot[:]

                # ignore qc-fails. this has to happen after the above, or the
                # previous agi's histogram won't be entered
                if not passed_qc:
                    print read_name
                    continue
                norm_hit_plot.append(100.0 * float(pos) / float(length - 20))
                last_agi_id = agi_id
            these_hits = hits_cur.fetchmany()
        cur.executemany(self.insert_histogram, histograms)
        self.sqlite_con.commit()
        cur.executemany(self.update_target_hits, updated_targets)
        self.sqlite_con.commit()
        cur.close()

    def run(self, sam_files):
        """
        Accepts a list of tuples:
            (sam_file_name, sample_type)
            sam_file_name: string which describes path to sam file
            sample_type: strings such as "WT", "Mutant" to group sample
        """
        for (sam_file_name, sample_type) in sam_files:
            sample_id = self.add_sample(sam_file_name, sample_type)
            print "Sample in %s (type %s) has id %i" %\
             (sam_file_name, sample_type, sample_id)
            self.fill_targets_table(sample_id)
            print "filled target table for sample %i" % sample_id
            self.fill_hits_table(sample_id)
            print "filled hits table for sample %i" % sample_id
            self.fill_histograms_table(sample_id)
            print "filled histogram table for sample %i" % sample_id
#            self.make_overall_histogram(sample_id)
#            print "made overall histogram for sample %i" % sample_id
#            self.make_culmulative_histogram(sample_id)
#            print "made culmulative histogram for sample %i" % sample_id
#            self.make_individual_histograms(sample_id)
#            print "made individual histograms for sample %i" % sample_id
        self.sqlite_con.close()
