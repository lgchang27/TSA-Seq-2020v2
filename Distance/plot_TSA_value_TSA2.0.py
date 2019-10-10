#!/usr/bin/python
# Programmer : Liguo Zhang, Yang Zhang
# Date: 
# Last-modified: 19 Feb 2018 07:41:00 PM

import os,sys,argparse
import numpy as np
import tempfile
from TSA_utility import *
import matplotlib
import matplotlib.pyplot as plt

def ParseArg():
    ''' This function parse the argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-w','--wig',type=str,dest="wig",help="wig file (e.g. smoothed TSA-seq signal wig file)")
    p.add_argument('-g','--gap',type=str,dest="gap",help="gap region file")
    return p.parse_args()

def LoadGap(fin_name):
    '''Load genome gaps'''
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

def CollectValue(wig_file,gap):
    """Collecct all 20kb-binnned values"""
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
    return array

def Main():
    global args
    args=ParseArg()
    gap_table = LoadGap(args.gap)
    value=CollectValue(args.wig,gap_table)
    max=np.max(value)
    min=np.min(value)
    print "max is: "+str(max)
    print "min is: "+str(min)

if __name__=="__main__":
    Main()
