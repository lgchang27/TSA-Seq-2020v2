#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 20180415
# Last-modified: 20190929

import os,sys,argparse
from TSA_utility import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-b','--bed',type=str,dest="bed",help="name of bed file (cell2 > cell1) of the changed domains")
    p.add_argument('-b2','--bed2',type=str,dest="bed2",help="name of bed file 2 (cell1 > cell2) of the changed domains")
    p.add_argument('-g1','--gf1',type=str,dest="gf1",help="name of gene expression file 1")
    p.add_argument('-g2','--gf2',type=str,dest="gf2",help="name of gene expression file 2")
    p.add_argument('-o','--output',type=str,dest="output",help="output file name")
    p.add_argument('-y','--yy',type=str,dest="yy",help="y axis indicator")
    return p.parse_args()

def getGene(gf,bedChrom,bedHashh):
    '''This function reads the file for gene expression (FPKM) results, identify all protein-coding genes within the changed domains from input, and return an array for all the genes containing FPKM value, region mean residual and gene ID. '''
    global args
    args=ParseArg()
    FPKMarray=[]
    for line in ReadFromFile(gf):
        row = line.strip().split()
        if row[0]=="gene_id":
            continue
        chrom=row[17]
        for chr in bedChrom:
            if chr==chrom:
                start=row[18]
                end=row[19]
                for region in bedHashh[chrom]:
                    if int(start)>=int(region[0]) and int(end)<=int(region[1]) and 'protein_coding' in row: # the whole gene must be within a domain to be called, and only report protein-coding genes
                        FPKM=row[6]
                        geneID=row[0]
                        FPKMarray.append((float(FPKM),float(region[2]),str(geneID)))
    return FPKMarray

def readBed(bedfile):
    '''Read bed file, identify all protein coding genes within the domains, and return a list of all FPKM log2-fold changes in the two cell lines'''
    bedChr=[]
    bedHash={}
    array=[]
    chrom=None
    for line in ReadFromFile(bedfile):
        row = line.strip().split()
        if row[0]!=chrom: 
            if len(array)!=0:
                bedHash[chrom]=array
                array=[]
            chrom = row[0]
            bedChr.append(chrom)
            array.append((row[1],row[2],row[3])) # domain start pos, domain end pos, domain mean TSA residual
        else:
            array.append((row[1],row[2],row[3])) # domain start pos, domain end pos, domain mean TSA residual
    bedHash[chrom]=array
    print bedChr
    chrom=None

    '''Get all protein-coding genes within the domains and their FPKM values from the two cell lines, region mean residuals, and gene ID'''
    FPKMarray1=getGene(args.gf1,bedChr,bedHash)
    FPKMarray2=getGene(args.gf2,bedChr,bedHash)
    print 'FPKMarray1: n='+str(len(FPKMarray1))
    print 'FPKMarray2: n='+str(len(FPKMarray2))
    if len(FPKMarray1) != len(FPKMarray2):
        print "The gene number in the two arrays are different"
        sys.exit("terminated")     
    
    '''Calculate FPKM log2-fold changes and return a list'''
    FPKM_fold_change=[]
    aa=0
    bb=0
    cc=0
    for i in range (len(FPKMarray1)):
        fpkmFoldChange=np.log2((FPKMarray1[i][0]+0.1)/(FPKMarray2[i][0]+0.1))
        FPKM_fold_change.append(fpkmFoldChange)  
    return FPKM_fold_change 

def Main():
    global args
    args=ParseArg()
    fpkm1=readBed (args.bed)
    fpkm2=readBed (args.bed2)
    
    '''kernel density plot of FPKM log2-fold changes'''
    plt.subplots(figsize=(2,8))
    ax = sns.kdeplot(fpkm1,shade=True,color="blue",vertical=True)
    ax = sns.kdeplot(fpkm2,shade=True,color="red",vertical=True)
    #plt.ylim (-15,15) # adjust range
    plt.xlabel (args.yy,multialignment='center',fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=15, width=1.5)
    plt.subplots_adjust(bottom=.2, left=.5)
    plt.savefig(args.output+"density.eps",format='eps')

if __name__=="__main__":
    Main()
