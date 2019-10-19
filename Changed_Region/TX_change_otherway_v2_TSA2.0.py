#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 20180415
# Last-modified: 20190307

import os,sys,argparse
from TSA_utility import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-b','--bed',type=str,dest="bed",help="name of bed file of the changed domains")
    p.add_argument('-g1','--gf1',type=str,dest="gf1",help="name of gene expression file 1")
    p.add_argument('-g2','--gf2',type=str,dest="gf2",help="name of gene expression file 2")
    p.add_argument('-o','--output',type=str,dest="output",help="output file name")
    p.add_argument('-y','--yy',type=str,dest="yy",help="y axis indicator")
    p.add_argument('-geneID','--geneID',type=str,dest="geneID",help="output file name for gene id")
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
                        FPKMarray.append((float(FPKM),float(region[2]),str(geneID))) # return FPKM value, region mean residual, and gene ID
    return FPKMarray

def Main():
    global args
    args=ParseArg()
    geneIDout=WriteToFile(args.geneID)
    
    '''Read bed file, and return a dictionary with chromosomes as keys and the changed domains of each chromosomes'''
    bedChr=[]
    bedHash={}
    array=[]
    chrom=None
    for line in ReadFromFile(args.bed):
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
    
    '''Calculate FPKM log2-fold changes and collect the rescaled TSA-Seq score residuals'''
    FPKM_fold_change=[]
    percentile_change=[]
    aa=0
    bb=0
    for i in range (len(FPKMarray1)):
        fpkmFoldChange=np.log2((FPKMarray1[i][0]+0.1)/(FPKMarray2[i][0]+0.1))
        FPKM_fold_change.append(fpkmFoldChange)
        if FPKMarray1[i][1]==FPKMarray2[i][1]:
            percentileChange=FPKMarray1[i][1]
            percentile_change.append(percentileChange*(-1))  #the cell line comparison to identify chnaged domains is directional, so here change the negative residuals to positive values for plotting
        if fpkmFoldChange < 0:
            if FPKMarray1[i][2]!=FPKMarray2[i][2]:
                print "The gene id in the two arrays is different"
            print >>geneIDout, '%s' % (str(FPKMarray1[i][2]))
            aa=aa+1
        if fpkmFoldChange == 0:
            bb=bb+1
    
    '''Report gene numbers'''
    print 'biased expressed gene number: '+str(aa) 
    print 'unbiased expressed gene number: '+str(bb) 
    print 'FPKM_fold_change: n='+str(len(FPKM_fold_change))
    print 'percentile_change: n='+str(len(percentile_change))
    if len(FPKM_fold_change) != len(percentile_change):
        print "The numbers in the two arrays are different"
        sys.exit("terminated")   

    '''Scatter plot'''
    plt.scatter(percentile_change,FPKM_fold_change,marker=',', color='black', s=2)
    #plt.xlim (10,30)   #reset range
    #plt.ylim (-15,15)  #reset range  
    plt.xlabel('Domain mean TSA change\n'+args.output,fontsize=15)       
    plt.ylabel ('Protein-coding gene FPKM\nlog2 ('+args.yy+')',multialignment='center',fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=15, width=1.5)
    plt.subplot().spines['left'].set_linewidth(1.5)
    plt.subplot().spines['bottom'].set_linewidth(1.5)
    plt.subplot().spines['right'].set_linewidth(1.5)
    plt.subplot().spines['top'].set_linewidth(1.5)
    plt.subplots_adjust(bottom=.2, left=.2)
    plt.savefig(args.output+".eps",format='eps')


if __name__=="__main__":
    Main()
