#!/usr/bin/python
# Programmer : Liguo Zhang, modified from Yang Zhang's code
# Date: 
# Last-modified: 26 Feb 2018 11:38 PM

import os,sys,argparse
import numpy as np
import tempfile
from TSA_utility import *

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-w','--wig',type=str,dest="wig",help="wig file (e.g. smoothed TSA-seq signal wig file)")
    p.add_argument('-q','--quantile',type=int,dest="quantile",help="quantile segment number, 10 means divide 100%% into 10 segment")
    p.add_argument('--quantilelist',type=str,dest="quantilelist",help="quantile list, eg, 0,10,20,25,45,50,100, the first one must be 0 and the last one must be 100 and should be in an ascending order, the delimiter must be comma")
    p.add_argument('-g','--gap',type=str,dest="gap",help="gap region file")
    p.add_argument('-o','--output',type=str,dest="output",help="output file name without .x")
    p.add_argument('-gg','--genome',type=str,dest="genome",help="human genome file (two column, first col is chromosome name, second chromosome is chromosome size")
    return p.parse_args()

def LoadGap(fin_name):
    '''This function will load gap regions that will be removed for later normalization'''
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

def FindQuantile(wig_file,quantile_num,gap,is_list=False):
    """This function will find quantile segment upper and lower limits"""
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
        if num % 10000 == 0:
            logging("Process: %d line processed" % (num))
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
    if is_list is False:
        quantile_list = (range(0,100,100/quantile_num) + [100])[1:]
    else:
        quantile_list = sorted([int(item.strip()) for item in quantile_num.split(',')])
        if 0 not in quantile_list or 100 not in quantile_list:
            error("0 and 100 must exist in quantilelist")
            exit(1)
        quantile_list = quantile_list[1:]
    return CombineList(quantile_list, np.percentile(array,quantile_list))

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
    '''This function will read wig files and convert TSA-Seq scors to percentile and then write into wig and bedbraph files'''
    num = 1
    for line in ReadFromFile(wig_file):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        if row[0] == 'track':
            continue
        if num % 10000 == 0:
            logging("Process: %d line processed" % (num))
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
    if args.quantile is not None:
        quantile_seg = FindQuantile(args.wig,args.quantile,gap_table)
    elif args.quantilelist is not None:
        quantile_seg = FindQuantile(args.wig,args.quantilelist,gap_table,True)
    else:
        error("quantile or quantilelist must be set")
        exit(1)
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
