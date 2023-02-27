#!/usr/bin/env python3
import os
import sys
import os.path as path
import argparse
import math


def get_base_fast(genome_reference_file_handle, index_dict, chrom, start, end=0):
    '''
    Get the specific bases from an index genome file. (samtools faidx <ref.fasta>)
    For the speed we need to pass the genome_reference_file_handle and index dictionary

    Parameter:
        **genome_reference_file_handle**: file IO
            a genome reference fasta file open handle

        **index_dict**: dict
            index is dictionary in faidx format, can be read from an index genome file (samtools faidx <ref.fasta>)
            {'1': [249250621, 52, 60, 61],
            '2': [243199373, 253404903, 60, 61],
            '3': [198022430, 500657651, 60, 61], .... }

        **chrom**: str
            chromosome name

        **start**: int
            start position (1-based)

        **end**: optional, int
            end position (1-based)

    Return:
        Bases from chromosome
    '''

    start_offset_int = -1
    end_offset_int = -1

    # in-chrome offset + in-line offset + in-file offset
    if start % index_dict[chrom][2] != 0:
        start_offset_int = math.floor(start / index_dict[chrom][2]) * index_dict[chrom][3] + start % index_dict[chrom][2] + index_dict[chrom][1]
    else:
        start_offset_int = math.floor((start / index_dict[chrom][2] - 1)) * index_dict[chrom][3] + index_dict[chrom][2] + index_dict[chrom][1]

    genome_reference_file_handle.seek(start_offset_int - 1)

    if end == 0:
        end_offset_int = start_offset_int
    else:
        if end % index_dict[chrom][2] != 0:
            end_offset_int = math.floor(end / index_dict[chrom][2]) * index_dict[chrom][3] + end % index_dict[chrom][2] + index_dict[chrom][1]
        else:
            end_offset_int = math.floor((end / index_dict[chrom][2] - 1)) * index_dict[chrom][3] + index_dict[chrom][2] + index_dict[chrom][1]

    try:
        seq_str = genome_reference_file_handle.read(end_offset_int - start_offset_int + 1).replace('\n', '')
    except Exception as ex:
        print(ex)
        return('')

    return(seq_str)

    def __get_offset(index_dict, chrom, pos):
        # in-chrome offset + in-line offset + in-file offset
        if pos % index_dict[chrom][2] != 0:
            n = math.floor(pos / index_dict[chrom][2]) * index_dict[chrom][3] + pos % index_dict[chrom][2] + index_dict[chrom][1]
        else:
            n = math.floor((pos / index_dict[chrom][2] - 1)) * index_dict[chrom][3] + index_dict[chrom][2] + index_dict[chrom][1]
        return(n)

    offset_start_int = __get_offset(genome_file, chrom, start)
    if end != 0:
        offset_end_int = __get_offset(genome_file, chrom, end)
    else:
        offset_end_int = offset_start_int

    if offset_end_int < offset_start_int:
        message = '\n\tThe end position (now {0}) must be greater than start position (now {1}).'.format(end, start)
        raise IndexError(message)

    with open(file=path.realpath(path.expanduser(genome_file)), mode='rt') as f:
        f.seek(offset_start_int - 1)
        raw_str = f.read(offset_end_int - offset_start_int + 1)
        return(raw_str.replace('\n', ''))
