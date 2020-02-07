#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 
# Last-modified: 20191215, calculate mean and sd distance for each percentile or maxminNorm score
# wig1 is percentile/maxminNorm score wig file, wig2 is distance wig file

import os,sys,argparse
import numpy as np
import tempfile
from TSA_utility import *
import matplotlib
import matplotlib.pyplot as plt

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-w1','--wig1',type=str,dest="wig1",help="percentile wig file1 (e.g. wig file1)")
    p.add_argument('-w2','--wig2',type=str,dest="wig2",help="percentile wig file2 (e.g. wig file2)")    
    p.add_argument('-o','--out',type=str,dest="output",help="output name")
    return p.parse_args()

def CollectValue(wigfile):
    hash={}
    array=[]
    num=1
    all=[]
    for line in ReadFromFile(wigfile):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        if row[0] == 'track':
            continue
        if num % 10000 == 0:
            logging("Process: %d line processed in wig" % (num))
        num += 1
        if row[0] == 'variableStep': # new track chrom begin
            if len(array)!=0:
                if chrom != row[1].split('=')[1]:
                    hash[chrom]=array
                    array=[]
            chrom = row[1].split('=')[1]
        else:
            array.append((row[0],row[1]))
            all.append(float(row[1]))
        hash[chrom]=array
    print "total line num:"+str(num)
    print "all value num:"+str(len(all))
    print "max vallue:"+str(max(all))
    print "min vallue:"+str(min(all))
    return hash

def Main():
    global args
    args=ParseArg()
    TSApercentile1=CollectValue(args.wig1)
    TSApercentile2=CollectValue(args.wig2)
    chrs_1 = TSApercentile1.keys()
    chrs_2 = TSApercentile2.keys()
    chrs = list(set(chrs_1) & set(chrs_2))
    print chrs
    print 'chr num:'+str(len(chrs))
    TSApercentile1_array=[]
    TSApercentile2_array=[]
    for i in range(len(chrs)):
        A_list = TSApercentile1[chrs[i]]
        B_list = TSApercentile2[chrs[i]]       
        if len(A_list) != len(B_list):
            print "The bin number for chrs[i] in the two tsa files are different"
            sys.exit("terminated") 
        for j in range (len(A_list)):
            if A_list[j][0]==B_list[j][0]:
                TSApercentile1_array.append(float(A_list[j][1]))
                TSApercentile2_array.append(float(B_list[j][1]))
    print max(TSApercentile1_array)
    print min(TSApercentile1_array)
    print max(TSApercentile2_array)
    print min(TSApercentile2_array)
    print "TSApercentile1 bin num:"+str(len(TSApercentile1_array))
    print "TSApercentile2 bin num:"+str(len(TSApercentile2_array))
    
    dic={}
    out = WriteToFile(args.output)
    print >>out, "%s\t%s\t%s\t%s" % ('TSA-score','n', 'mean','std')
    for i in range (100):
        pp=i+1
        dic['perc'+str(pp)]=[]
        for j in range (len(TSApercentile1_array)):
            if int(TSApercentile1_array[j])+50 == int(pp):
                dic['perc'+str(pp)].append(TSApercentile2_array[j])
        n=len(dic['perc'+str(pp)])
        mean=np.mean(dic['perc'+str(pp)])
        std=np.std(dic['perc'+str(pp)])
        print >>out, "%d\t%d\t%.6f\t%.6f" % (pp,n,mean,std)
    
if __name__=="__main__":
    Main()
