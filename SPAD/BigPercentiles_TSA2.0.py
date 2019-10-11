#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 20180419
# Last-modified: 
# This code will identify bins with percentile above a certain threshold and then merge adjacent bins to call regions. This code will also generate histograms showing statistics of region size and number.

import os,sys,argparse
from TSA_utility import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-w','--wig',type=str,dest="wig",help="wig file (e.g. TSA-seq percentile wig file)")
    p.add_argument('-o','--output',type=str,dest="output",help="output tsa compare file name")
    p.add_argument('-p','--percentileThreshold', type=float, dest='percentileThreshold',help='set the threshold for identify regions above it')
    p.add_argument('-g','--genome',type=str,dest="genome",help="human genome file (two column, first col is chromosome name, second chromosome is chromosome size")
    p.add_argument('-win','--window',type=int,dest="window",help="window size of input wigs")
    return p.parse_args()

def identify(wig,hgenome):
    '''This function will identify all bins with percentile above a certain threshold'''
    global args
    hash={}
    array=[]
    num=1
    bedout = WriteToFile(args.output + '_above_'+str(args.percentileThreshold)+'.bed')
    for line in ReadFromFile(wig):
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
            array.append((row[0],float(row[1])+50))
    hash[chrom]=array
    print num
    chrs = hash.keys()
    print "chrs:"
    print chrs
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
    array=[]
    identifyHash={}
    for chrom in chrs:       
        chrom_size = int(LoadGenome(hgenome)[chrom])
        A_list = hash[chrom]
        for j in range (len(A_list)):            
                percentile=float(A_list[j][1])
                if percentile >= args.percentileThreshold:
                    if j== len(A_list) - 1:
                        print >>bedout, "%s\t%d\t%d\t%.6f" % (chrom,int(A_list[j][0]),chrom_size,float(percentile))
                    else:
                        print >>bedout, "%s\t%d\t%d\t%.6f" % (chrom,int(A_list[j][0]),int(A_list[j][0])+args.window,float(percentile))
                    array.append((A_list[j][0],percentile))
        identifyHash[chrom]=array
        array=[]
    bedout.flush() 
    return identifyHash

def MergeSegment_adjacent (filename,outname):
    """This function will merge adjacent bins to call a region"""
    out = WriteToFile(outname)
    chrom = None
    sizeList=[]
    percList=[]
    for line in ReadFromFile(filename):
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        if row[0] != chrom:
            if chrom is not None:
                print >>out, "%s\t%s\t%s\t%.6f" % (chrom,start,stop,np.mean(percentileList))
                sizeList.append(int(stop)-int(start))
                percList.append(abs(np.mean(percentileList)))
            chrom = row[0]
            start = row[1]
            stop = row[2]
            percentileList=[float(row[3])]
        else:
            if row[1] != stop:
                print >>out, "%s\t%s\t%s\t%.6f" % (chrom,start,stop,np.mean(percentileList))
                sizeList.append(int(stop)-int(start))
                percList.append(abs(np.mean(percentileList)))
                chrom = row[0]
                start = row[1]
                stop = row[2]
                percentileList=[float(row[3])]
            else:
                stop = row[2]
                percentileList.append(float(row[3]))
    print >>out, "%s\t%s\t%s\t%.6f" % (chrom,start,stop,np.mean(percentileList))
    sizeList.append(int(stop)-int(start))
    percList.append(abs(np.mean(percentileList)))
    out.flush()
    return sizeList,percList

def Main():
    global args
    args=ParseArg()
    genome = LoadGenome(args.genome)
    aboveThreshold=identify(args.wig, args.genome)
    chrs = aboveThreshold.keys()
    print chrs
    print 'chr num:'+str(len(chrs))
    array=[]
    for i in range(len(chrs)):
        list = aboveThreshold[chrs[i]]      
        for j in range (len(list)):
                array.append(list[j][1])
    print "TSA bin num above percentile threshold:"+str(len(array))
    print "max TSA percentile above threshold:"+str(max(array))
    print "min TSA percentile above threshold:"+str(min(array))
    plt.hist(array, bins=np.arange(np.floor(min(array))-5, np.ceil(max(array))+5, 1), facecolor='k')
    plt.xlabel('TSA_percentile_above_threshold')
    plt.ylabel ('count')
    plt.title('Histogram of bins with TSA percentile above threshold')
    plt.savefig(args.output+"_bin_above_threshold.png")
    plt.close()

    SortBed(args.output + '_above_'+str(args.percentileThreshold)+'.bed')
    mergename = args.output + '_above_'+str(args.percentileThreshold)+"_mergeAdjacent.bed"
    size,perc= MergeSegment_adjacent(args.output + '_above_'+str(args.percentileThreshold)+'.bed',mergename)

    sort="sort -k1,1 -k2,2n %s > %s" % (mergename, args.output + '_above_'+str(args.percentileThreshold)+"_mergeAdjacent.sorted.bed")
    bed2bb = "./utilities/bedToBigBed %s %s %s" % (args.output + '_above_'+str(args.percentileThreshold)+"_mergeAdjacent.sorted.bed", args.genome, args.output + '_above_'+str(args.percentileThreshold)+"_mergeAdjacent.bb")
    os.system(sort)
    os.system(bed2bb)

    plt.hist(size, bins=np.arange(np.floor(min(size))-20000, np.ceil(max(size))+20000, 20000), facecolor='k')
    plt.xlabel('TSA_percentile_above_threshold_domain_size')
    plt.ylabel ('count')
    plt.title('Histogram of TSA_percentile_above_threshold_domain_size')
    plt.annotate('domain number:'+str(len(size))+'\nmean domian size:'+str(int(np.mean(size)))+'\nmedian domian size:'+str(int(np.median(size)))+'\nmax domian size:'+str(int(np.max(size)))+'\nmin domian size:'+str(int(np.min(size))),xy=(0.5, 0.6), xycoords='axes fraction',fontsize=10)
    plt.savefig(args.output+'_above_'+str(args.percentileThreshold)+"_domain_size_hst.png")
    plt.close()

    plt.hist(perc, bins=np.arange(np.floor(min(perc))-5, np.ceil(max(perc))+5, 1), facecolor='k')
    plt.xlabel('TSA_percentile_above_threshold_domain')
    plt.ylabel ('count')
    plt.title('Histogram of TSA_percentile_above_threshold_domain')
    plt.annotate('mean perc:'+str(int(np.mean(perc)))+'\nmedian perc:'+str(int(np.median(perc)))+'\nmax perc:'+str(int(np.max(perc)))+'\nmin perc:'+str(int(np.min(perc))), xy=(0.5, 0.6), xycoords='axes fraction',fontsize=10)
    plt.savefig(args.output+"_domain_hst.png")
    plt.close()    

    logging("DONE!!!")

if __name__=="__main__":
    Main()
