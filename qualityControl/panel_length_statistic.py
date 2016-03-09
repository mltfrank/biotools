#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
    Statistic script for gene panel length.
    Author: Wang Bingchen
    Date: 2016/3/7
"""

import os
import sys
import logging

format_list = ['table', 'csv']

def usage():
    print """
    Analysis a interval or bed file to statistic the distribution of gene panel's length.
    Print the length of a gene panel, following with a percentage of gene panel whose length is longer than this value.
    Usage:
        python panel_length_statistic.py -I <input_file> [-L <low_bound>] [-F <format>] [-Q <length>]
            -I <input_file>: Bed file or interval file, writen in text.
            -L <low_bound>: Optional. Integer value. length with more percentage than 'low_bound' will be print. [0.8]
            -F <format>: Optional. Can be 'table' or 'csv' ['table']
            -Q <length>: Optional. When -Q is set, this program print the percentage of gene panel with length longer than query length. Query length must be a positive integer
          """

def parse_arguments():
    """
        Parse arguments to get settings
    """
    config = {}
    try:
        config['format'] = 'table'
        config['low_bound'] = 0.8
        config['query'] = -1
        for idx, arg in enumerate(sys.argv):
            if arg == "-I":
                arg_val = sys.argv[idx+1]
                config['input_file'] = arg_val
            elif arg == '-L':
                arg_val = sys.argv[idx+1]
                config['low_bound'] = float(arg_val)
            elif arg == '-F':
                arg_val = sys.argv[idx+1]
                config['format'] = arg_val
            elif arg == '-Q':
                arg_val = sys.argv[idx+1]
                config['query'] = int(arg_val)
        check_flag = True
        if not config.get('input_file', None):
            logging.error("Must have a '-I' option")
            check_flag = False
        if not config['format'] in format_list:
            logging.error("Unknown output format")
            check_flag = False
        if not check_flag:
            usage()
            sys.exit(2)
    except Exception, e:
        logging.error("Arguments parse error")
        print e
        sys.exit(2)
    return config

def statistic(input_file_path):
    """
        Traverse all lines in interval or bed file.
        Split each line by '\t' and treat colume 2 and 3 as the start and end coordinate of panel.
        Put the result into a map.
        Arguments:
            input_file_path: interval or bed file path
        Return:
            Gene panel count in this file.
            A list from length to gene panel count, sorted by length decrease.
    """
    count_map = {}
    panel_count = 0
    with open(input_file_path) as input_file:
        for line in input_file:
            if line.startswith("@"):
                continue
            try:
                line_split = line.split('\t')
                length = int(line_split[2]) - int(line_split[1])
                if length < 0:
                    length = 0-length
                cur_count = count_map.get(length, 0)
                count_map[length] = cur_count+1
                panel_count += 1
            except:
                # just jump over line
                logging.info("Jump over line '%s'"%(line[:-1]))
                continue
    return (panel_count, sorted(count_map.iteritems(), key=lambda d:d[0], reverse = True))

def cal_percentage(panel_count, count_list):
    """
        Calculate the percentage for each length.
        Arguments:
            count_list: list of tuple (length, count), with length decrease.
            panel_count: all count of panel.
        Return:
            a list from length to percentage.
    """
    percentage_list = []
    total_count = 0
    for (length, count) in count_list:
        total_count += count
        percentage = float(total_count) / panel_count
        percentage_list.append((length, percentage))
    return percentage_list

def print_result(percentage_list, low_bound, format, query):
    """
        Print result.
        Arguments:
            percentage_list: list with tuple (length, percentage)
            low_bound: program will print all length with percentage less than the bound.
            format: output format
    """
    inter_ch = '\t'
    if format == 'table':
        inter_ch = '\t'
    if format == 'csv':
        inter_ch = ','
    if query >= 0:
        last_length = 0
        last_percentage = 0
        print "length%spercentage"%(inter_ch)
        for (length, percentage) in percentage_list:
            if length >= query:
                last_length = length
                last_percentage = percentage
            else:
                print "%d%s%.10f"%(last_length, inter_ch, last_percentage)
                return
        print inter_ch.join((str(query), '1'))
    else:
        for (length, percentage) in percentage_list:
            print inter_ch.join((str(length), str(percentage)))
            if percentage > low_bound:
                break


def main():
    config = parse_arguments()
    (panel_count, count_list) = statistic(config['input_file'])
    percentage_list = cal_percentage(panel_count, count_list)
    print_result(percentage_list, config['low_bound'], config['format'], config['query'])

if __name__ == "__main__":
    logging.basicConfig(format='%(asctime)s -- %(levelname)s:%(message)s', level = logging.INFO)
    main()
