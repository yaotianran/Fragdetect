#!/usr/bin/env python3
# python scripts/get_GC.py genome.reference
import sys
import time
import re
import os.path as path
from utils import utils


def main(argvList = sys.argv, argv_int = len(sys.argv)):
    '''
    # input temp/pos.tsv
    # chr1	10399	10598	+
    # chr1	10464	10516	+

    # the output is a one-column headless data frame , each row is a fragment
    # output temp/gc.tsv
    # column one is GC ratio (float),
    '''
    t = time.time()
    INPUT = path.realpath(path.expanduser('temp/pos.tsv'))
    OUTPUT = path.realpath(path.expanduser('temp/gc.tsv'))
    refer = path.realpath(path.expanduser(argvList[1]))

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

    genome_reference_file_handle = open(refer, 'rt')
    prog = re.compile('G|C')
    gc_float = 0.0
    with open(INPUT, mode = 'rt') as in_f, open(OUTPUT, mode = 'wt', buffering = 100000) as out_f:
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
                seq_str = utils.get_base_fast(genome_reference_file_handle, index_dict, chrom, start, end)
                gc_float = len(prog.findall(seq_str)) / (end - start + 1)
            except:
                gc_float = 'NA'
            out_f.writelines(str(gc_float) + '\n')

    genome_reference_file_handle.close()

    print('{i}: {time:.2f} seconds done'.format(i = i, time = time.time() - t))
    return


if __name__ == '__main__':
    # refer = '/home/user/sda1/data/human_g1k_v37/human_g1k_v37.fasta'
    # r = __get_gc_ratio(genome_file = refer, chrom = '3', start = 100000, end = 100100)

    main()

