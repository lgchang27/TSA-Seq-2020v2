#!/usr/bin/python
# Programmer : Liguo Zhang, Yang Zhang
# Date: 
# Last update:20181024

import os,sys,argparse,math
import numpy as np
import tempfile
from TSA_utility import *

def parse_argument():
    ''' This Function Parse the Argument '''
    p = argparse.ArgumentParser()
    p.add_argument('-i','--input',dest = "input",type = str, help = "the name of input file")
    p.add_argument('-o','--output',dest = "output",type = str, help = "the name of output file")
    p.add_argument('-g','--genome',dest = "genome",type = str, help = "the name of genome file")
    p.add_argument('-gap','--gap',dest = "gap",type = str, help = "the name of gap file")
    p.add_argument('-y0','--y0',dest = "y0",type = float, help = "y0 value")
    p.add_argument('-A','--A',dest = "A",type = float, help = "A value")
    p.add_argument('-R0','--R0',dest = "R0",type = float, help = "R0 value")   	
    return p.parse_args()

global args
args = parse_argument()

def LoadGap(fin_name):
    '''load genome gap regions'''
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


def Main():
    wout = WriteToFile(args.output + '.wig')
    gap = LoadGap(args.gap)
    chrom = None
    span = None
    start = None
    num = 1
    for line in ReadFromFile(args.input):
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
            print >>wout, "variableStep chrom=%s span=%d" % (chrom,span)
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
                overlap_len += min((start+span), gap[chrom][mm][1]) - max(gap[chrom][mm][0], start) # multiple gap regions insides each bin
                #if float( min((start + span),gap[chrom][mm][1]) - max(gap[chrom][mm][0], start) )/float(span) < 0.75:
                #    continue
                #is_overlap = True
            assert overlap_len <= span
            if float(overlap_len)/float(span) >= 0.75:
                is_overlap = True
            if not is_overlap:
                y = math.pow(2,float(value))
                if y-args.y0 >= 1e-6 and y < (float(args.y0)+float (args.A)):
                    x = 1/(args.R0)*math.log((y-(args.y0))/(args.A)) 
                elif y-args.y0 < 1e-6:
                    x = 1/(args.R0)*math.log((1e-6)/(args.A))
                else:
                    x = 0.0
                print >>wout, "%d\t%.6f" % (start,x)
    print "total line processed:"+str(num-1)
    wout.flush()
    wig2bw = "utilities/wigToBigWig -clip %s %s %s" % ( args.output+".wig" , args.genome, args.output + '.bw')
    os.system(wig2bw)
    logging("DONE!!!")

if __name__=="__main__":
    Main()              		