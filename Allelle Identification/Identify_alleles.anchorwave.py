#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author  : Qionghou Li
# @FileName: Identify_alleles.anchorwave.py

##version 2.3 use seqlite3 to make the process short

import sys
import sqlite3


# create table dat(qstart int, qend int, sstart int, send int, qlens int, slens int,identity int, qstand text , sstand text , qchr text, schr text);
# .separator "\t"
# .import HXSA_vs_HXSB.gmap.f.maf.coords dat
#
#
# create table dat(qstart INTEGER, qend INTEGER, sstart INTEGER, send INTEGER, qlens INTEGER, slens INTEGER,identity INTEGER, qstand text , sstand text , qchr text, schr text);
# .separator "\t"
# .import HXSA_vs_HXSB.gmap.f.maf.coords dat
# select * from dat where qchr="Chr01A" and schr="Chr01B" and qstart>24956814 and qend<24967253 ;
#
# create table dat(g1 text, g2 text, simarity INTEGER, length INTEGER, mismatch INTEGER, gapopen INTEGER,qstart INTEGER, qend INTEGER , sstart INTEGER , send INTEGER, evalue INTEGER, bitscore INTEGER);
# .separator "\t"
# .import HXSA_vs_HXSB.ultrasens.e1.nomax.blast dat


class ReadBed(object):
    def __init__(self, fn):
        self.fn = fn
        dct = self._read()
        new_dct = self._read2(fn)
        self.dct = dct
        self.new_dct = new_dct

    def _read(self):
        fn = self.fn
        dct = {}
        with open(fn, 'r') as op1:
            op1.seek(0, 0)
            for i in op1:
                dct[i.strip('\n').split('\t')[3]] = i
        return dct

    def get(self, ids):
        dct = self.dct
        if ids in dct.keys():
            return dct[ids]
        else:
            sys.exit('gene ID was not in your bed file'.format(ids))

    def get_full(self, ids):
        dct = self.dct
        if ids in dct.keys():
            tmp = dct[ids].strip('\n').split('\t')
            return tmp[0], tmp[1], tmp[2]
        else:
            print(ids)

    def _read2(self, fn):
        dct = {}
        with open(fn, 'r') as op1:
            op1.seek(0, 0)
            for i in op1:
                ln = i.strip('\n').split('\t')
                genes = gene(ln[3], ln[0], ln[1], ln[2], ln[5])
                dct[ln[3]] = genes
        return dct


class gene(object):
    def __init__(self, ids, chr, s, e, stand):
        self.chr = chr
        self.ids = ids
        self.s = int(s)
        self.e = int(e)
        self.stand = stand
        self.full = [chr, int(s), int(e), ids, stand]


##gene class will given,loc1=[s,e+1],loc2=[s,e+1]
def loc_index(g1, bed, v):
    # print(g1, v)
    closest, id = '', ''
    loc1 = [int(bed.new_dct[g1].s), int(bed.new_dct[g1].e)]
    for i in v:
        ## judge overlap
        loc2 = [int(bed.new_dct[i].s), int(bed.new_dct[i].e)]
        mean1, mean2 = map(lambda x: np.mean(x), [loc1, loc2])
        if closest == '':
            closest = abs(mean1 - mean2)
            id = i
        elif abs(mean1 - mean2) <= closest:
            closest = abs(mean1 - mean2)
            id = i
    return g1, id, closest


def judge_similarity(blastdb, gene1, gene2, Maxthreshold=100, Minthreshold=0.2):
    blastdb.row_factory = dict_factory
    blasts = blastdb.cursor()
    blasts.execute('select * from dat where g1=\"{}\" and g2=\"{}\";'.format(gene1, gene2))
    lst = blasts.fetchall()
    if len(lst) > 1:
        print(gene1, gene2, lst)
    elif len(lst) == 1:
        if lst[0]['mismatch'] == 0 and lst[0]['gapopen'] == 0 and lst[0]['simarity'] == int(Maxthreshold):
            return 1
        else:
            return 0
    else:
        print('{} and {} pairs not show in diamond'.format(gene1, gene2))


def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d


def judge_region(con, qchr, qstart, qend, tchr, tstart, tend):
    con.row_factory = dict_factory
    cur = con.cursor()
    qstart, qend, tstart, tend = int(qstart), int(qend), int(tstart), int(tend)

    # cur.execute(
    #     'select * from dat where qchr=\"{}\" and schr=\"{}\" and qstart<={} and qend>={} and sstart <= {} and send >={};'.format(
    #         qchr, tchr, int(qstart), int(qend), int(tstart), int(tend)))
    cur.execute(
        'select * from dat where qchr=\"{}\" and schr=\"{}\" and (qstart<={} and qend>={} OR qstart>={} and qstart<={} OR qend>={} and qend<={}) and (sstart<={} and send>={} OR sstart>={} and sstart<={} OR send>={} and send<={});'.format(
            qchr, tchr, int(qstart), int(qend), qstart, qend, qstart, qend, tstart, tend, tstart, tend, int(tstart),
            int(tend)))
    lst = cur.fetchall()
    return lst


if __name__ == '__main__':
    anchors, anchorwave, bed, blast = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    fbed = ReadBed(bed)
    blastdb = sqlite3.connect(blast)
    con = sqlite3.connect(anchorwave)
    dct1 = {}
    #  dct2={}
    with open(anchors, 'r') as op1, open(sys.argv[5], 'w') as op2, open('NotinAnchorwave.genePairs.txt', 'w') as op3:
        for i in op1:
            if i.startswith('#'):
                continue
            else:
                lne = i.strip('\n').split('\t')
                dct1.setdefault(lne[0], []).append(lne[1])
        ##直接判断是否为allele，如果是就留下，如果都是就都留下
        for k, v in dct1.items():
            qchr, qstart, qend = fbed.new_dct[k].chr, fbed.new_dct[k].s, fbed.new_dct[k].e
            #   print(qchr,qstart,qend)
            for i in v:
                tchr, tstart, tend = fbed.new_dct[i].chr, fbed.new_dct[i].s, fbed.new_dct[i].e
                #      print(qchr,qchr[:-1],tchr)
                #  print(tchr,tstart,tend)
                if qchr[:-1] == tchr[:-1]:
                    lst = judge_region(con, qchr, qstart, qend, tchr, tstart, tend)
                    #  print(lst)
                    if len(lst) > 0:
                        if judge_similarity(blastdb, k, i):
                            op2.write(k + '\t' + i +'\t'+'1'+ '\n')
                        else:
                            op2.write(k + '\t' + i + '\n')
                    else:
                        print("{}\t{}".format(k, i), file=op3)
                        print('{} and {} not show in anchorwave file'.format(k, i))
                    #    dct2[k]=i
                ## add blast filter steps
                else:
                    continue
