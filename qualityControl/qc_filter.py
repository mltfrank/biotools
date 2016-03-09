#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
    QC filter for fastq.
    A filted fastq file and a report file will be generated after execution.
    Author: Wang Bingchen
    Date: 2016/3/7
"""

import os
import sys
import logging

from numpy import average

def usage():
    print """
    Usage:
    A script to generate qc filted fastq file and filt report.
    Strategy to filt read is described as follow:
    We got a read with 4 lines, include 2 description line, a base line and a quality line.
    We traverse the quality line from tail to head, if we found a series of base with high quality, \
we keep the sub string from head to the series of base as high quality clipped read. Then we do \
different filter on this read. A LOW_AVERAGE_QUAL filter and SHORT_READ_FILTER is used currently.
    We record the base-pair and read count we filted in report.

    python qc_filter.py -I <input_fastq> -O <output_fastq> [-R <report_file>] [-S <score_type>] [-TQ <tail_qual>] [-NT <num_of_high_qual_in_tail>] [-ML <min_read_length>] [-AQ <average_quality>]
      -I <input_fastq>: required, input fastq
      -O <output_fastq>: required, filted output fastq
      -R <report_file>: optional, output report file. If not set, no report file will be generate
      -S <score_type>: optional, quality score type in input fastq file, can be 'sanger', 'solexa', 'illumina1.3+', 'illumina1.5+', 'illumina1.8+' ['sanger']
      -TQ <tail_qual>: optional, threshold of high tail quality [30]
      -NT <num_of_high_qual_in_tail>: if count of high quality tail base is more than num_of_high_qual_in_tail, cut read there.[10]
      -ML <min_read_length>: min length of a clipped read for filter. Set according to the length distribution of selected gene panel. [50]
      -AQ <average_quality>: average quality threshold of filter. [30]
    """


class LowAverageQualFilter(object):
    @classmethod
    def filt(cls, read, config):
        """
            Filter to see if the average quality of each base is lower than min_average_qual.
            Arguments:
                read: qual_array, quality for each base pair.
                config: min_average_qual, threshold.
            Return:
                If this read can be filted by this filter.
        """
        qual_array = [ord(qual) for qual in read[3].strip()]
        return average(qual_array) >= add_quality_offset(config['min_average_qual'], config['score_type'])

class ShortReadFilter(object):
    @classmethod
    def filt(cls, read, config):
        """
            Filter to see if a read is shorter than min_read_length.
            Arguments:
                read: base_array, base array
                config: min_read_length, threshold
            Return:
                if this read can be filted by this filter
        """
        # minus one to skip '\n'
        return len(read[1])-1 >= config['min_read_length']

filter_list = [ShortReadFilter, LowAverageQualFilter]
def filt_read(read, config):
    """
        Filt using filter in filter_list
        If add a new filter in the futhure, just add a filter class and add into the filter_list.
        Arguments:
            read: 4 element list, for each line in a read.
            config: config dict.
        Return:
            (if_pass, filter_id)
    """
    for filter_idx, filter_cls in enumerate(filter_list):
        if not filter_cls.filt(read, config):
            return (False, filter_idx)
    return (True, -1)

def clip_low_qual_tail(read, tail_qual, num):
    """
        Read clipper to cut the low quality tail of read.
        If met num base-pair with quality higher than tail_qual, stop cutting.
        Arguments:
            read: a list with 4 element list, each list represents desc_array, base_array, desc_array2 and qual_array.
            tail_qual: tail quality threshold.
            num: number of expected high quality base
        Return:
            clipped read length
    """
    if len(read) != 4:
        logging.error("Illegal read")
        sys.exit(2)
    qual_array = [ord(qual) for qual in read[3].strip()]
    length = len(qual_array)
    high_count = 0
    cut_idx = length-1
    for i in range(length):
        idx_reverse = length -i -1
        if qual_array[idx_reverse] >= tail_qual:
            high_count += 1
            if high_count > num:
                break
        else:
            high_count = 0
            cut_idx = idx_reverse-1
    return cut_idx+1

def check_read_format(read):
    """
        Check read format.
        For speed, just check the start charactor of two discription line.
        Return if read format is right.
    """
    try:
        return read[0][0] == '@' and read[2][0] == '+'
    except:
        return False

def clip_filt_reads(input_fastq, output_fastq, report_file, config):
    """
        Clip low quality tail.
        Filt reads.
        Generate report if needed.
        Arguments:
            input_fastq : input fastq
            output_fastq: output fastq
            report_file: report file path
            config: read filter and clipper configuration in this dict.
    """
    total_clipped_base = 0
    total_kept_base = 0
    total_filted_read = [0 for cls in filter_list]
    total_kept_read = 0
    line_count = 0

    generate_report = config['generate_report']
    fixed_tail_qual = add_quality_offset(config['tail_qual'], config['score_type'])
    num_of_high_qual_in_tail = config['num_of_high_qual_in_tail']

    with open(input_fastq, 'r') as input_file:
        with open(output_fastq, 'w') as output_file:
            next_line = input_file.readline()
            while(next_line):
                # read a read from file
                read = []
                read.append(next_line)
                read.append(input_file.readline())
                read.append(input_file.readline())
                read.append(input_file.readline())
                if not check_read_format(read):
                    logging.error("Illegal read format in line %d"%(line_count+1))
                    sys.exit(2)
                line_count += 4

                # clip low quality read end.
                raw_length = len(read[1]) -1
                clipped_length = clip_low_qual_tail(read, fixed_tail_qual, num_of_high_qual_in_tail)
                read[1] = read[1][:clipped_length]+'\n'
                read[3] = read[3][:clipped_length]+'\n'

                # filt low quality read
                (if_pass, filt_idx) = filt_read(read, config)
                #if if_pass:
                if if_pass:
                    output_file.write(read[0])
                    output_file.write(read[1])
                    output_file.write(read[2])
                    output_file.write(read[3])

                # generate report statistic
                if generate_report:
                    total_clipped_base += raw_length - clipped_length
                    total_kept_base += clipped_length
                    if not if_pass:
                        total_filted_read[filt_idx] += 1
                    else:
                        total_kept_read += 1

                next_line = input_file.readline()

    # mark a dict to generate report
    if generate_report:
        generate_report_file(report_file, total_clipped_base, total_kept_base,
                filter_list, total_filted_read, total_kept_read)

def generate_report_file(report_file_path, total_clipped_base, total_kept_base,
        filter_list, total_filted_read, total_kept_read):
    """
        Generate report file
        Arguments:
            total_filted_read: is a list for each filter.
    """
    total_base = total_clipped_base + total_kept_base
    clip_base_rate = float(total_clipped_base) / float(total_base)
    kept_base_rate = float(total_kept_base) / float(total_base)
    clip_kept_base_rate = float(total_clipped_base) / float(total_kept_base)

    filted_read_sum = sum(total_filted_read)
    total_read = filted_read_sum + total_kept_read
    filt_read_rate = float(filted_read_sum) / float(total_read)
    kept_read_rate = float(total_kept_read) / float(total_read)
    filt_kept_read_rate = float(filted_read_sum) / float(total_kept_read)
    with open(report_file_path, 'w') as report_file:
        report_file.write("Total clipped base: %d\t\tclip rate: %.4f%%\n"%(total_clipped_base, clip_base_rate*100))
        report_file.write("Total kept base: %d\t\tkept rate: %.4f%%\n"%(total_kept_base, kept_base_rate*100))
        report_file.write("Clipped base / Kept base: %.4f%%\n"%(clip_kept_base_rate*100))
        report_file.write("\n")
        report_file.write("Total filted read: %d\t\tfilt rate: %.4f%%\n"%(filted_read_sum, filt_read_rate*100))
        report_file.write("Total kept read: %d\t\tkept rate: %.4f%%\n"%(total_kept_read, kept_read_rate*100))
        report_file.write("Filted read / Kept read: %.4f%%\n"%(filt_kept_read_rate*100))
        for cls, filted_read in zip(filter_list, total_filted_read):
            filter_rate = float(filted_read) / float(total_read)
            report_file.write("\tFilter '%s': %d\t\tfilt rate:%.4f%%\n"%(cls.__name__, filted_read, filter_rate*100))

def add_quality_offset(old_threshold, score_type):
    """
        Define offset for different fastq score type.
        Arguments:
            old_threshold: old_threshold
        Return:
            Fixed quality threshold.
    """
    if score_type == 'sanger':
        return old_threshold + 33
    if score_type == 'solexa' or score_type == 'illumina1.3+' \
            or score_type == 'illumina1.5+' or score_type == 'illumina1.8+':
        return old_threshold + 64
    else:
        logging.error("Unsupported score type for fastq.")
        sys.exit(2)

def parse_arguments():
    config = {}
    config['generate_report'] = False
    config['score_type'] = 'sanger'
    config['tail_qual'] = 30
    config['num_of_high_qual_in_tail'] = 10
    config['min_read_length'] = 50
    config['min_average_qual'] = 30
    try:
        for idx, arg in enumerate(sys.argv):
            if arg == '-I':
                config['input_fastq'] = sys.argv[idx+1]
            elif arg == '-O':
                config['output_fastq'] = sys.argv[idx+1]
            elif arg == '-R':
                config['report_file'] = sys.argv[idx+1]
                config['generate_report'] = True
            elif arg == '-S':
                config['score_type'] = sys.argv[idx+1]
            elif arg == '-TQ':
                config['tail_qual'] = int(sys.argv[idx+1])
            elif arg == '-NT':
                config['num_of_high_qual_in_tail'] = int(sys.argv[idx+1])
            elif arg == '-ML':
                config['min_read_length'] = int(sys.argv[idx+1])
            elif arg == '-AQ':
                config['min_average_qual'] = int(sys.argv[idx+1])
        check_arg = True
        if not config.get('input_fastq', None):
            logging.error("Must provide a input fastq")
            check_arg = False
        if not config.get('output_fastq', None):
            logging.error("Must provide a output fastq")
            check_arg = False
        if not check_arg:
            usage()
            sys.exit(2)
    except Exception, e:
        logging.error("Error when parse arguments")
        print e
        usage()
        sys.exit(2)
    return config

def main():
    config = parse_arguments()
    clip_filt_reads(config['input_fastq'], config['output_fastq'], config['report_file'], config)


if __name__ == '__main__':
    main()
