#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 20190307
# Last-modified: 

'''This code will compare two cell lines with two replicates each, take input mean and standard diviation from fitted gaussian distribution, and return p-values for each 20kb bins '''

import os,sys,argparse
from TSA_utility import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from TSA_utility import *
import scipy.stats
from scipy.stats import norm
from math import log

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-c1r1','--c1r1',type=str,dest="c1r1",help="wig file1 (e.g. TSA-seq percentile wig file1, cell line1 replicate1)")
    p.add_argument('-c1r2','--c1r2',type=str,dest="c1r2",help="wig file1 (e.g. TSA-seq percentile wig file1, cell line1 replicate2)")
    p.add_argument('-c2r1','--c2r1',type=str,dest="c2r1",help="wig file1 (e.g. TSA-seq percentile wig file1, cell line2 replicate1)")
    p.add_argument('-c2r2','--c2r2',type=str,dest="c2r2",help="wig file1 (e.g. TSA-seq percentile wig file1, cell line2 replicate2)")
    p.add_argument('-o','--output',type=str,dest="output",help="output tsa compare file name")
    p.add_argument('-g','--genome',type=str,dest="genome",help="human genome file (two column, first col is chromosome name, second chromosome is chromosome size")
    p.add_argument('-w','--window',type=int,dest="window",help="window size of input wigs")
    p.add_argument('-mean','--mean',type=float,dest="mean",help="mean of gaussian distribution")
    p.add_argument('-std','--std', type=float, dest='std',help='std of gaussian distribution')
    p.add_argument('-maxRsd','--maxRsd', type=float, dest='maxRsd',help='max of residue')
    p.add_argument('-minRsd','--minRsd', type=float, dest='minRsd',help='min of residue')
    p.add_argument('-t','--t', type=float, dest='threshold',help='threshold')
    return p.parse_args()

def CollectValue(wigfile):
    ''' This Function Parse the Argument '''
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
        num += 1
        if row[0] == 'variableStep': 
            if len(array)!=0:
                if chrom != row[1].split('=')[1]:
                    hash[chrom]=array
                    array=[]
            chrom = row[1].split('=')[1]
        else:
            array.append((row[0],row[1]))
            all.append(float(row[1]))
        hash[chrom]=array
    print "all value num:"+str(len(all))
    print "max vallue:"+str(max(all))
    print "min vallue:"+str(min(all))
    return hash

def Compare(first_wig,second_wig, third_wig,fourth_wig, hgenome, win):
    '''This function will take input mean and standard diviation from fitted gaussian distribution, and return p-values for each 20kb bins, output a wig file'''
    global args
    first_hash=CollectValue(first_wig)
    second_hash=CollectValue(second_wig)
    third_hash=CollectValue(third_wig)
    fourth_hash=CollectValue(fourth_wig)
    wigout = WriteToFile(args.output + '.wig')
    chrs_A = first_hash.keys()
    chrs_B = second_hash.keys()
    chrs_C = third_hash.keys()
    chrs_D = fourth_hash.keys()
    if len(chrs_A) != len(chrs_C):
        print "The chromsome number between two cell lines are different, are you sure continue? (Yes/No)"
        go_on = raw_input("> ")
        if go_on == "No":
            sys.exit("interupt by user")
    chrs = list(set(chrs_A) & set(chrs_C))
    sort_table = {}
    for chrom in chrs:
        if chrom[0:3].lower() == 'chr':
            try:
                sort_table[int(chrom[3:])] = chrom
            except ValueError:
                sort_table[chrom[3:]] = chrom
        else:
            sort_table[chrom] = chrom
    sorted_chrom = sorted(sort_table.keys())
    newchrs=[]
    for name in sorted_chrom:
        newchrs.append(sort_table[name]) 
    chrs=newchrs
    print "sorted chrs:"
    print chrs
    for chrom in chrs:       
        chrom_size = int(LoadGenome(hgenome)[chrom])
        A_list = first_hash[chrom]
        B_list = second_hash[chrom] 
        C_list = third_hash[chrom]
        D_list = fourth_hash[chrom]       
        print >>wigout, "variableStep chrom=%s span=%d" % (chrom,win)
        if not ((len(A_list) == len(B_list) and len(A_list) == len(C_list) and len(A_list) == len(D_list))):
            print "The bin number for chrs[i] 2 tsa files are different"
            sys.exit("terminated")  
        for j in range (len(A_list)):
            if (A_list[j][0]==B_list[j][0] and A_list[j][0]==C_list[j][0] and A_list[j][0]==D_list[j][0]):
                change=0.5*((float(C_list[j][1])+float(D_list[j][1]))-(float(A_list[j][1])+float(B_list[j][1])))
                if change >= 0:
                    P=scipy.stats.norm(args.mean,args.std).cdf(-change)
                if change < 0:
                    P=scipy.stats.norm(args.mean,args.std).cdf(change)             
                if j== len(A_list) - 1 and (chrom_size-int(A_list[j][0])) < win:
                    print >>wigout, "variableStep chrom=%s span=%d" % (chrom,chrom_size-int(A_list[j][0]))
                    print >>wigout, "%d\t%.6f" % (int(A_list[j][0]),float((-log(P,10))))
                else:
                    print >>wigout, "%d\t%.6f" % (int(A_list[j][0]),float((-log(P,10))))
    wigout.flush()
    wig2bw = "utilities/wigToBigWig -clip %s %s %s" % ( args.output+'.wig', args.genome, args.output+'.bw')
    os.system(wig2bw) 

def Main():
    global args
    args=ParseArg()
    genome = LoadGenome(args.genome)
    
    ###test
    Pthrshd=1-scipy.stats.norm(args.mean,args.std).cdf(args.threshold)
    print "test: p value for threshold is:"+str(Pthrshd)
    Pmean=1-scipy.stats.norm(args.mean,args.std).cdf(args.mean)
    print "test: p value for mean is:"+str(Pmean)
    PmaxRsd=scipy.stats.norm(args.mean,args.std).cdf(-args.maxRsd)
    print "test: p value for maxRsd is:"+str(PmaxRsd)
    PminRsd=scipy.stats.norm(args.mean,args.std).cdf(args.minRsd)
    print "test: p value for minRsd is:"+str(PminRsd)

    Compare(args.c1r1,args.c1r2, args.c2r1,args.c2r2, args.genome,args.window)

    logging("DONE!!!")

if __name__=="__main__":
    Main()
