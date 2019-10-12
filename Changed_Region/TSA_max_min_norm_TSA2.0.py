#!/usr/bin/python
# Programmer : Liguo Zhang, modified form Yang's code
# Date: 
# Last-modified: 26 Feb 2018 11:38 PM

import os,sys,argparse
import numpy as np
import tempfile
from TSA_utility import *

'''This code will normalize TSA-Seq enrichment scores based on max and min values'''
'''If set quantile_num=100, this normalizaton is equal to the funcion: 
    Scaled enrichment score (bin i) = (TSA-Seq enrichment score (bin i) - min) / (max - min) * 100
    (min assigned to 1 instead of 0)'''
'''The term "quantile" in this code doesn't mean real quantiles of values. Instead, it is just used as a term.'''


def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-w','--wig',type=str,dest="wig",help="wig file (e.g. smoothed TSA-seq signal wig file)")
    p.add_argument('-q','--quantile',type=int,dest="quantile",help="quantile segment number, 10 means divide 100%% into 10 segment")
    p.add_argument('-g','--gap',type=str,dest="gap",help="gap region file")
    p.add_argument('-o','--output',type=str,dest="output",help="output file name without .x")
    p.add_argument('-max','--max',type=float,dest='max', help='normalized max tsa value')
    p.add_argument('-min','--min',type=float,dest='min', help='normalized min tsa value')
    p.add_argument('-gg','--genome',type=str,dest="genome",help="human genome file (two column, first col is chromosome name, second chromosome is chromosome size")
    return p.parse_args()

def LoadGap(fin_name):
    '''This function will load genome gap regions that will be removed from the normalization'''
    table = {}
    for line in ReadFromFile(args.gap):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        chrom = row[0]
        start = int(row[1])+1 # convert to 1-base
        stop = int(row[2])
        if chrom not in table.keys():
            table[chrom] = []
        table[chrom].append((start,stop))
    return table

def FindQuantile(wig_file,quantile_num,gap, TSAmax,TSAmin):
    '''This function will find upper and lower limits for segments to assign TSA-Seq enrichment scores for normalization'''
    '''if quantile_num=100, this normalizaton is equal to the funcion: 
    Scaled enrichment score (bin i) = (TSA-Seq enrichment score (bin i) - min) / (max - min) * 100
    (min assigned to 1 instead of 0)'''
    chrom = None
    span = None
    start = None
    array = []
    num = 1
    for line in ReadFromFile(wig_file):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        if row[0] == 'track':
            continue
        num += 1
        if row[0] == 'variableStep': # new track chrom begin
            chrom = row[1].split('=')[1]
            span = int(row[2].split('=')[1])
            start = None
        else:
            is_overlap = False
            try:
                start = int(row[0])
                value = float(row[1])
            except ValueError:
                warning("Warning: Read wig error. Unknown line format for line:%s" % (line.strip()))
            overlap_len = 0
            for mm in range(len(gap[chrom])):
                if start + span <= gap[chrom][mm][0]:
                    continue
                if start > gap[chrom][mm][1]:
                    continue
                overlap_len += min((start+span), gap[chrom][mm][1]) - max(gap[chrom][mm][0], start) 
            assert overlap_len <= span
            if float(overlap_len)/float(span) >= 0.75:
                is_overlap = True
            if not is_overlap:
                array.append(value)
    quantile_list = (range(0,100,100/quantile_num) + [100])[1:]
    cutofflist=np.arange(TSAmin,TSAmax,(TSAmax-(TSAmin))/quantile_num)[1:] 
    TSArmax=np.max(array)
    print np.max(array)
    print np.min(array)
    if len (cutofflist)==len (quantile_list):
        cutofflist=cutofflist[:-1]
    cutofflist=np.append(cutofflist,TSArmax+0.1)
    return CombineList(quantile_list, cutofflist)

def CombineList(list_a,list_b):
    result = []
    if len(list_a) != len(list_b):
        print len(list_a)
        print len(list_b)
        error("Can't combine two lists with different length")
        exit(1)
    for nn in range(len(list_a)):
        result.append((list_a[nn],list_b[nn]))
    return result

def WriteToSeg(out,wout,wig_file,gap,quantile):
    '''This function will assign each TSA-Seq enrichment scores (20kb bin) into corresponding groups for normalization'''
    num = 1
    for line in ReadFromFile(wig_file):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        if row[0] == 'track':
            continue
        num += 1
        if row[0] == 'variableStep':
            chrom = row[1].split('=')[1]
            span = int(row[2].split('=')[1])
            start = None
            print >>wout, "variableStep chrom=%s span=%d" % (chrom,span)
        else:
            is_overlap = False
            try:
                start = int(row[0])
                value = float(row[1])
            except ValueError:
                warning("Parsing wig file error: Unknown line format for line: %s" % (line.strip()))
            overlap_len = 0
            for mm in range(len(gap[chrom])):
                if start + span <= gap[chrom][mm][0]:
                    continue
                if start > gap[chrom][mm][1]:
                    continue
                overlap_len += min((start+span), gap[chrom][mm][1]) - max(gap[chrom][mm][0], start)
            assert overlap_len <= span
            if float(overlap_len)/float(span) >= 0.75:
                is_overlap = True
                print >>out, chrom + '\t' + str(start-1) + '\t' + str(start+span-1) + '\t' + 'N' + '\t' + '{:.6f}'.format(value)
            if is_overlap:
                continue
            else:
                for nn in range(len(quantile)):
                    if value < quantile[nn][1]:
                        break
                print >>wout, "%d\t%.1f" % (start,quantile[nn][0]-50)
                print >>out, chrom + '\t' + str(start-1) + '\t' + str(start+span-1) + '\t' + str(quantile[nn][0]) + '\t' + '{:.6f}'.format(value)

def Main():
    global args
    args=ParseArg()
    out = WriteToFile(args.output+'.bedgraph')
    outwig = WriteToFile(args.output + '.wig')
    gap_table = LoadGap(args.gap)
    quantile_seg = FindQuantile(args.wig,args.quantile,gap_table,args.max,args.min)
    logging("Results: Dividing upper bound is " + '\t'.join('{:.3f}'.format(item[1]) for item in quantile_seg))
    WriteToSeg(out,outwig,args.wig,gap_table,quantile_seg)
    out.flush()
    outwig.flush()
    SortBed(args.output+'.bedgraph')
    wig2bw = "utilities/wigToBigWig -clip %s %s %s" % ( args.output+".wig" , args.genome, args.output + '.bw')
    os.system(wig2bw)
    logging("Finish: TSA_quantile DONE!!!")

if __name__=="__main__":
    Main()
