#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 20180415
# Last-modified: 20191218

import os,sys,argparse
from TSA_utility import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-b01','--bed01',type=str,dest="bed01",help="bed file 01")
    p.add_argument('-b02','--bed02',type=str,dest="bed02",help="bed file 02")
    p.add_argument('-b03','--bed03',type=str,dest="bed03",help="bed file 03")
    p.add_argument('-g1','--gtf1',type=str,dest="gtf1",help="gtf file1")
    p.add_argument('-g2','--gtf2',type=str,dest="gtf2",help="gtf file2")
    p.add_argument('-g3','--gtf3',type=str,dest="gtf3",help="gtf file3")
    p.add_argument('-g0','--gtf0',type=str,dest="gtf0",help="gtf file0")
    p.add_argument('-o','--output',type=str,dest="output",help="output file name")
    return p.parse_args()

def getGene(gtf,bedChrom,bedHashh):
    '''This function reads the file for gene expression (FPKM) results, identify all protein-coding genes within the changed domains from input, and return an array for all the genes containing FPKM value and gene ID. '''
    global args
    args=ParseArg()
    FPKMarray=[]
    for line in ReadFromFile(gtf):
        row = line.strip().split()
        if row[0]=="gene_id":
            continue
        chrom=row[17]
        for chr in bedChrom:
            if chr==chrom:
                start=row[18]
                end=row[19]
                for region in bedHashh[chrom]:
                    if int(start)>=int(region[0]) and int(end)<=int(region[1]) and 'protein_coding' in row:
                        FPKM=row[6]
                        geneID=row[0]
                        FPKMarray.append((float(FPKM),str(geneID)))
    return FPKMarray

def readBed (bedf):
    '''read bed file'''
    bedChr=[]
    bedHash={}
    array=[]
    chrom=None
    for line in ReadFromFile(bedf):
        row = line.strip().split()
        if row[0]!=chrom:  #new chrom begain
            if len(array)!=0:
                bedHash[chrom]=array
                array=[]
            chrom = row[0]
            bedChr.append(chrom)
            array.append((row[1],row[2]))  # domain start pos, domain end pos
        else:
            array.append((row[1],row[2]))  # domain start pos, domain end pos
    bedHash[chrom]=array
    return bedHash,bedChr

def Main():
    global args
    args=ParseArg()
    regions01,chrs01 = readBed(args.bed01)
    regions02,chrs02 = readBed(args.bed02)
    regions03,chrs03 = readBed(args.bed03)

    FPKMarray01a=getGene(args.gtf0,chrs01,regions01)
    FPKMarray01b=getGene(args.gtf1,chrs01,regions01)
    if len(FPKMarray01a) != len(FPKMarray01b):
        print "The gene number in the two arrays are different"
        sys.exit("terminated")     
    FPKM_fold_change01=[]
    aa=0
    bb=0
    for i in range (len(FPKMarray01a)):
        fpkmFoldChange=np.log2((FPKMarray01a[i][0]+0.1)/(FPKMarray01b[i][0]+0.1))
        FPKM_fold_change01.append(fpkmFoldChange)
        if FPKMarray01a[i][1]!=FPKMarray01b[i][1]:
            print "The gene in the two arrays are different"
            sys.exit("terminated")     
        if fpkmFoldChange > 0:
            aa=aa+1
        if fpkmFoldChange == 0:
            bb=bb+1
    print 'biased expressed gene number of bed01: '+str(aa) 
    print 'unbiased expressed gene number of bed01: '+str(bb) 
    print 'FPKM_fold_change: n='+str(len(FPKM_fold_change01))



    FPKMarray02a=getGene(args.gtf0,chrs02,regions02)
    FPKMarray02b=getGene(args.gtf2,chrs02,regions02)
    if len(FPKMarray02a) != len(FPKMarray02b):
        print "The gene number in the two arrays are different"
        sys.exit("terminated")     
    FPKM_fold_change02=[]
    aa=0
    bb=0
    for i in range (len(FPKMarray02a)):
        fpkmFoldChange=np.log2((FPKMarray02a[i][0]+0.1)/(FPKMarray02b[i][0]+0.1))
        FPKM_fold_change02.append(fpkmFoldChange)
        if FPKMarray02a[i][1]!=FPKMarray02b[i][1]:
            print "The gene in the two arrays are different"
            sys.exit("terminated")     
        if fpkmFoldChange > 0:
            aa=aa+1
        if fpkmFoldChange == 0:
            bb=bb+1
    print 'biased expressed gene number of bed02: '+str(aa) 
    print 'unbiased expressed gene number of bed02: '+str(bb) 
    print 'FPKM_fold_change: n='+str(len(FPKM_fold_change02))



    FPKMarray03a=getGene(args.gtf0,chrs03,regions03)
    FPKMarray03b=getGene(args.gtf3,chrs03,regions03)
    if len(FPKMarray03a) != len(FPKMarray03b):
        print "The gene number in the two arrays are different"
        sys.exit("terminated")     
    FPKM_fold_change03=[]
    aa=0
    bb=0
    for i in range (len(FPKMarray03a)):
        fpkmFoldChange=np.log2((FPKMarray03a[i][0]+0.1)/(FPKMarray03b[i][0]+0.1))
        FPKM_fold_change03.append(fpkmFoldChange)
        if FPKMarray03a[i][1]!=FPKMarray03b[i][1]:
            print "The gene in the two arrays are different"
            sys.exit("terminated")     
        if fpkmFoldChange > 0:
            aa=aa+1
        if fpkmFoldChange == 0:
            bb=bb+1
    print 'biased expressed gene number of bed03: '+str(aa) 
    print 'unbiased expressed gene number of bed03: '+str(bb) 
    print 'FPKM_fold_change: n='+str(len(FPKM_fold_change03))

    sns.kdeplot(FPKM_fold_change01,shade=True,color="blue",vertical=True)
    sns.kdeplot(FPKM_fold_change02,shade=True,color="orange",vertical=True)
    sns.kdeplot(FPKM_fold_change03,shade=True,color="green",vertical=True)
    plt.ylim (-12,12)    #adjust range
    plt.tick_params(axis='both', which='major', labelsize=15, width=1.5)
    plt.subplot().spines['left'].set_linewidth(1.5)
    plt.subplot().spines['bottom'].set_linewidth(1.5)
    plt.subplot().spines['right'].set_linewidth(1.5)
    plt.subplot().spines['top'].set_linewidth(1.5)    
    plt.subplots_adjust(bottom=.2, left=.5)
    plt.savefig(args.output+"density.eps",format='eps')
    plt.close()

if __name__=="__main__":
    Main()
