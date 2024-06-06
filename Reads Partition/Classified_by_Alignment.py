#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author  : Qionghou Li
# @FileName: Classified_by_Alignment.py

import sys
import logging
import pysam
from numpy import mean

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s')


def check_sorted(bam, sorted='queryname'):
    header = dict(bam.header)
    if 'HD' in header:
        if sorted == 'name':
            if header['HD']['SO'] == 'queryname':
                print('sorted OK')
                return 1
        else:
            if header['HD']['SO'] == 'coordinate':
                raise ValueError("Please sort your bam with -n options, you are now coordinate")
                # sys.exit(-1)
    else:
        raise ValueError('Please sort your bam with -n options')
        # sys.exit(-1)


class ReadBam(object):
    def __init__(self, f, sorted='name'):
        # bamf = pysam.AlignmentFile(f, "rb")
        # self._check_sorted(bamf, sorted)
        self.bam = f

    def get_paired_scores(self):
        bam = self.bam
        dct = {}
        for i in bam:
            reads, scores = i.query_name, i.get_tag('AS')
            dct.setdefault(reads, []).append(scores)
        return dct


class CompareBams(object):
    # bam1 and bam2 are two dicts generated from ReadBam.get_paired_scores()
    def __init__(self, bam1, bam2):
        set1, set2 = set(bam1.keys()), set(bam2.keys())
        read1 = set1 - set2
        read2 = set2 - set1
        read_1_2 = set1 & set2
        self.read1 = read1
        self.read2 = read2
        self.read_1_2 = read_1_2
        self.bam1 = bam1
        self.bam2 = bam2

    def _judge_reads(self, bam1, bam2, i):
        v1, v2 = bam1.get(i), bam2.get(i)
        # print(v1,v2)
        mean1, mean2 = mean(list(map(int, v1))), mean(list(map(int, v2)))
        if mean1 > mean2:
            return 'a'
        if mean1 < mean2:
            return 'b'
        if mean1 == mean2:
            return 'c'

    def classified(self):
        read_1_2 = self.read_1_2
        op1 = open('{}.hap1.id'.format(sys.argv[1]), 'w')
        op2 = open('{}.hap2.id'.format(sys.argv[2]), 'w')
        hap1_cont, hap2_count, hap1_2 = 0, 0, 0
        bam1, bam2 = self.bam1, self.bam2
        for i in read_1_2:
            values = self._judge_reads(bam1, bam2, i)
            if values == 'a':
                hap1_cont += 1
                op1.write(i + '\n')
                # lst1.append(i)
                # print(i)
            elif values == 'b':
                hap2_count += 1
                op2.write(i + '\n')
                # lst2.append(i)
            # print(i)
            elif values == 'c':
                hap1_2 += 1
                op1.write(i + '\n')
                op2.write(i + '\n')
                # lst1.append(i)
                # lst2.append(i)
            # print(i)
        return hap1_cont, hap2_count, hap1_2


if __name__ == '__main__':
    bam1 = pysam.AlignmentFile(sys.argv[1], "rb")
    bam2 = pysam.AlignmentFile(sys.argv[2], "rb")
    # check_sorted(bam1)
    # check_sorted(bam2)
    # print('finish check')
    logging.debug(
        "Now Finish Check format steps"
    )
    dct1 = ReadBam(bam1).get_paired_scores()
    dct2 = ReadBam(bam2).get_paired_scores()
    # print('finish classified')
    logging.debug(
        "Now Finish Bam reading"
    )
    hap1_count, hap2_count, hap1_2 = CompareBams(dct1, dct2).classified()
    hap1_all = hap1_count + hap1_2
    hap2_all = hap2_count + hap1_2

    ratio1 = (hap1_all / hap2_all)
    ratio2 = hap1_count / hap2_count


    print(
        "\nHap1 reads {}\nHap2 reads {}\ncommon reads {}\nHap1 specifc reads {}\nHap2 specifc reads {}\n hap1_vs_hap2_ratio {:.2f}\nspecific_ratio_hap1_vs_hap2 {:.2f}\n".format(
            hap1_all, hap2_all, hap1_2, hap1_count, hap2_count, ratio1, ratio2)
    )
