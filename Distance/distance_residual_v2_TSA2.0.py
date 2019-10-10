#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 20181023
# Last-modified: 
# This code will compare two distance wig files, retun distance residual wig and bw file and a histgram for absolute residuals

import os,sys,argparse
from TSA_utility import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from TSA_utility import *

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-w1','--wig1',type=str,dest="wig1",help="wig file1 (e.g. calibrated wig file1)")
    p.add_argument('-w2','--wig2',type=str,dest="wig2",help="wig file2 (e.g. calibrated wig file2)")
    p.add_argument('-o','--output',type=str,dest="output",help="distance residual wig file as output")
    p.add_argument('-gap','--gap',type=str,dest="gap",help="gap region file")
    p.add_argument('-g','--genome',type=str,dest="genome",help="human genome file (two column, first col is chromosome name, second chromosome is chromosome size")
    p.add_argument('-w','--window',type=int,dest="window",help="window size of input wigs")
    return p.parse_args()

def LoadGap(fin_name):
    ''' Load genome gap regions'''
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

def Compare(first_wig,second_wig, hgenome, win,gap):
    '''This function take two wig files and calculate residuals for all values (20kb bin)'''
    
    global args
    first_hash={}
    second_hash={}
    
    '''collect values from first wig file'''
    array=[]
    num=1
    wigout = WriteToFile(args.output + '.wig')
    for line in ReadFromFile(first_wig):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        if row[0] == 'track':
            continue
        if num % 10000 == 0:
            logging("Process: %d line processed in wig1" % (num))
        num += 1
        if row[0] == 'variableStep': # new track chrom begin
            if len(array)!=0:
                if chrom != row[1].split('=')[1]:
                    first_hash[chrom]=array
                    array=[]
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
                array.append((row[0],row[1]))            
    first_hash[chrom]=array
    print "1st wigtotal line number:"+str(num-1)
    
    '''collect values from second wig file'''
    array=[]
    num=1
    for line in ReadFromFile(second_wig):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        if row[0] == 'track':
            continue
        if num % 10000 == 0:
            logging("Process: %d line processed in wig2" % (num))
        num += 1
        if row[0] == 'variableStep': # new track chrom begin
            if len(array)!=0:
                if chrom != row[1].split('=')[1]:
                    second_hash[chrom]=array
                    array=[]
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
                array.append((row[0],row[1]))            
    second_hash[chrom]=array
    print "2nd wig total line number:"+str(num-1)

    '''calculate absolute residuals between values from the two wig files'''
    chrs_A = first_hash.keys()
    chrs_B = second_hash.keys()
    if len(chrs_A) != len(chrs_B):
        print "The chromsome number between two wig files are different, are you sure continue? (Yes/No)"
        go_on = raw_input("> ")
        if go_on == "No":
            sys.exit("interupt by user")
    chrs = list(set(chrs_A) & set(chrs_B))
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
    compare_hash = {}
    array=[]
    for chrom in chrs:       
        chrom_size = int(LoadGenome(hgenome)[chrom])
        A_list = first_hash[chrom]
        B_list = second_hash[chrom]       
        print >>wigout, "variableStep chrom=%s span=%d" % (chrom,win)
        if len(A_list) != len(B_list):
            print "The bin number for chrs[i] in the two tsa files are different"
            sys.exit("terminated") 
        for j in range (len(A_list)):
            if A_list[j][0]==B_list[j][0]:
                change=abs(float(B_list[j][1])-float(A_list[j][1]))
                if j== len(A_list) - 1:
                    if chrom_size-int(A_list[j][0]) < args.window:
                        print >>wigout, "variableStep chrom=%s span=%d" % (chrom,chrom_size-int(A_list[j][0]))
                        print >>wigout, "%d\t%.6f" % (int(A_list[j][0]),float(change))
                    else:
                        print >>wigout, "%d\t%.6f" % (int(A_list[j][0]),float(change))
                else:
                    print >>wigout, "%d\t%.6f" % (int(A_list[j][0]),float(change))            
                array.append((A_list[j][0],change))
        compare_hash[chrom]=array
        array=[]
    wigout.flush()
    wig2bw = "utilities/wigToBigWig -clip %s %s %s" % ( args.output+'.wig', args.genome, args.output+'.bw')
    os.system(wig2bw)   
    return compare_hash

def Main():
    global args
    args=ParseArg()
    genome = LoadGenome(args.genome)
    gap_table = LoadGap(args.gap)
    TSAcompare=Compare(args.wig1,args.wig2, args.genome,args.window, gap_table)
    chrs = TSAcompare.keys()
    print chrs
    print 'chr num:'+str(len(chrs))
    TSA_array=[]

    for i in range(len(chrs)):
        TSA_list = TSAcompare[chrs[i]]      
        for j in range (len(TSA_list)):
                TSA_array.append(TSA_list[j][1])                
    mean=np.mean(TSA_array)
    median=np.median(TSA_array)
    sd=np.std(TSA_array)
    max=np.max(TSA_array)
    min=np.min(TSA_array)
    print "TSA bin num:"+str(len(TSA_array))
    print "max distance residual:"+str(max)
    print "min distance residual:"+str(min)
    print "mean distance residual:"+str(mean)
    print "median distance residual:"+str(median)
    print "sd distance residual:"+str(sd)

    plt.hist(TSA_array, bins=np.arange(0, 0.35, 0.005), facecolor='k')
    plt.xlabel('distance residual (um)', fontsize=20)
    plt.ylabel ('count',fontsize=20)
    plt.title('Histogram of Abs Residual',fontsize=20)
    plt.annotate('mean:'+str(np.around(mean,decimals=3))+'\nmedian:'+str(np.around(median,decimals=3))+'\nsd:'+str(np.around(sd,decimals=3)),xy=(0.5, 0.6), xycoords='axes fraction',fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20, width=2)
    plt.tick_params(axis='both', which='minor', labelsize=20, width=2)
    plt.subplot().set_xticks([0,0.1,0.2,0.3])
    plt.subplot().spines['right'].set_visible(False)
    plt.subplot().spines['top'].set_visible(False)
    plt.subplot().spines['left'].set_linewidth(2)
    plt.subplot().spines['bottom'].set_linewidth(2)    
    plt.subplots_adjust(bottom=.2, left=.2)
    plt.savefig(args.output+"_hst.eps", format='eps')
    plt.close()

    logging("DONE!!!")

if __name__=="__main__":
    Main()
