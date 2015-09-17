"""
Remove all mouse-derived reads from a Human PDX .bam

Usage:
    mouse_filter.py -b BAMFILE -o OUT_FASTQ_STUB [-m MOUSE_VARS]

Options:
    -b BAMFILE           Human .bam alignment file
    -o OUT_FASTQ_STUB    Human fastq output file stub (e.g. 2015-123_human)
    -m MOUSE_VARS        Strain-specific variants
"""
import pysam
from docopt import docopt
import logging
import sys
import gzip
from bitwise_flags import flags
__author__ = 'jgrundst'


class BamParser:

    def __init__(self, bamfile=None, outfile=None, mouse_vars=None):
        self.handler = None
        self.logger = logging.getLogger(__name__)
        self.setup_logging()
        self.bam = None
        self.load_bam(bamfile)
        self.headers = None

        self.logger.info(
            "Input params - \nbamfile: {}\noutfiles: {}\nmouse_vars: {}".format(
                bamfile, outfile + '_(1/2).fq.gz', mouse_vars
            )
        )

    def setup_logging(self):
        if not self.logger.handlers:
            self.logger.setLevel(logging.DEBUG)
            self.handler = logging.FileHandler('MouseFilter.log')
            self.handler.setLevel(logging.DEBUG)
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - [%(levelname)s]: %(message)s')
            self.handler.setFormatter(formatter)
            self.logger.addHandler(self.handler)
            self.logger.info("MouseFilter object created")

    def load_bam(self, bamfile):
        try:
            self.bam = pysam.AlignmentFile(bamfile, "rb")
        except IOError:
            self.logger.error("Cannot find bamfile: {}".format(bamfile))
            self.logger.info("Exiting")
            sys.exit(1)

    def show_headers(self):
        self.headers = self.bam.header
        for record_type, records in self.headers.items():
            print record_type
            for i, record in enumerate(records):
                print "\t{},".format(i+1)
                if type(record) == str:
                    print "\t\t{}".format(record)
                elif type(record) == dict:
                    for field, value in record.items():
                        print "\t\t{}\t{}".format(field, value)


    def evaluate_pair(self, read1, read2):
        if read1.cigarstring and read2.cigarstring:
            if ((len(read1.cigar) == 1 and read1.cigar[0][0] == 0) and
                    (len(read2.cigar) == 1 and
                        read2.cigar[0][0] == 0)):
                pass
            else:
                self.buffer_reads(read1, read2)
        else:
            self.buffer_reads(read1, read2)

    def extract_human(self):
        """
        Traverse .bam file for primary alignment pairs, skipping over
        secondary alignments.  pass primary alignments to evaluator and
        fill the fastq buffers.  dump to fq.gz files
        :return:
        """
        self.logger.info("Extracting human reads to outfiles...")
        self.read1_buffer = ''
        self.read2_buffer = ''
        self.buffer_read_count = 0
        pair_count = 0
        # the first read is always the primary alignment
        read1 = self.bam.next()
        while read1:

            read2 = self.bam.next()
            # make sure we're looking at primary alignment of read2
            while ((bitwise_flag_check(read2, 'second_in_pair') is False) and
                   (bitwise_flag_check(read2, 'not_primary_alignment')) is True):
                read2 = self.bam.next()
            pair_count += 1
            if read1.query_name != read2.query_name:
                out_message = '''Pair #{}
                Read1 {}
                Read2 {}
                Don't have the same names.  Quitting...'''
                raise ValueError(out_message.format(pair_count, read1, read2))

            # decide if we keep or toss!
            self.evaluate_pair(read1, read2)

            # if we have self.buffer_size reads in the buffer, output to gzip
            if self.buffer_read_count == self.buffer_size:
                self.print_fq_buffers()

            # Move on to next pair. Is there another read1?
            try:
                read1 = self.bam.next()
                while ((bitwise_flag_check(read1, 'second_in_pair') is True) and
                       (bitwise_flag_check(read1, 'not_primary_alignment')) is True):
                    read1 = self.bam.next()
            except StopIteration:
                self.logger.info('No more reads, end of file')
                read1 = None
        self.fq1.close()
        self.fq2.close()
        self.logger.info("Counted {} sequence pairs".format(pair_count))
        self.logger.info("Final buffer dump of {} pairs".format(
            self.buffer_read_count))
        self.print_fq_buffers()

    def count_perfect_matches(self):
        c = 0
        for read in self.bam:
            if len(read.cigar) == 1 and read.cigar[0][0] == 0:
                c += 1
        print "{} perfectly aligned reads".format(c)


def bitwise_flag_check(read, flag_string):
    if read.flag & flags[flag_string] == flags[flag_string]:
        return True
    else:
        return False


def main():
    args = docopt(__doc__)
    print args
    bp = BamParser(bamfile=args['-b'], outfile=args['-o'])
    bp.extract_human()

if __name__ == '__main__':
    main()
