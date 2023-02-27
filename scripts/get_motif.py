#!/usr/bin/env python3
# python scripts/get_motif.py genome.reference upstream downstream
import sys
import time
import math
import os.path as path
import argparse
import collections
from utils import utils


def main(argvList = sys.argv, argv_int = len(sys.argv)):
    '''
    # input temp/pos.tsv
    # chr1	10399	10598	+
    # chr1	10464	10516	+

    # the output is a three-column headless data frame (X1, X2, X3) , each row is a motif statistics
    # output temp/motif.tsv
    # column one is motif sequence (char),
    # column two is motif count at 5-end (integer)
    # column three is motif count at 3-end (integer)
    '''
    t = time.time()

    INPUT = path.realpath(path.expanduser('temp/pos.tsv'))
    OUTPUT = path.realpath(path.expanduser('temp/motif.tsv'))

    refer = path.realpath(path.expanduser(argvList[1]))
    upstream = int(argvList[2])
    downstream = int(argvList[3])

    indexfile = refer + '.fai'
    index_dict = {}
    with open(file = indexfile, mode = 'rt') as index_f:
        for i, line in enumerate(index_f):
            line_lst = line.strip().split()
            try:
                index_dict[line_lst[0]] = list(map(int, line_lst[1:]))
            except Exception:
                message = 'File {0} is corrupted at line {1}. Position query may raise an IndexError.'.format(indexfile, i + 1)
                print(message)
                next

    motif_dict = collections.defaultdict(lambda: [0, 0])
    genome_reference_file_handle = open(refer, 'rt')
    with open(INPUT, mode = 'rt') as in_f:
        for i, line_str in enumerate(in_f):
            i += 1
            if i % 100000 == 0:
                print(i)

            line_lst = line_str.strip().split()
            chrom = line_lst[0]
            start = int(line_lst[1])
            end = int(line_lst[2])
            strand = line_lst[3]

            try:
                if strand == '+':
                    start_5_int = start - upstream
                    end_5_int = start + downstream - 1
                    start_3_int = end - downstream + 1
                    end_3_int = end + upstream
                else:
                    start_5_int = end - downstream + 1
                    end_5_int = end + upstream
                    start_3_int = start - upstream
                    end_3_int = start + downstream - 1

                motif_5_end_str = utils.get_base_fast(genome_reference_file_handle, index_dict, chrom, start_5_int, end_5_int)
                motif_3_end_str = utils.get_base_fast(genome_reference_file_handle, index_dict, chrom, start_3_int, end_3_int)

                motif_dict[motif_5_end_str][0] += 1
                motif_dict[motif_3_end_str][1] += 1
            except:
                continue

    genome_reference_file_handle.close()

    with open(OUTPUT, mode = 'wt') as out_f:
        for key_str, value_lst in motif_dict.items():
            line_str = '{}\t{}\t{}\n'.format(key_str, value_lst[0], value_lst[1])
            out_f.writelines(line_str)

    print('{i}: {time:.2f} seconds done'.format(i = i, time = time.time() - t))
    return


if __name__ == '__main__':
    # refer = '/home/user/sda1/data/human_g1k_v37/human_g1k_v37.fasta'
    # r = __get_gc_ratio(genome_file = refer, chrom = '3', start = 100000, end = 100100)

    main()

