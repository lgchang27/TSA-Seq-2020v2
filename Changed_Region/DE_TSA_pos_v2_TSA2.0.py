#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 20190430
# Last modified: 20190928 

import os,sys,argparse
from TSA_utility import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-b','--bed',type=str,dest="bed",help="bed file for identified changed domains")
    p.add_argument('-csv','--csv',type=str,dest="csv",help="csv file from DE analysis")
    p.add_argument('-o','--output',type=str,dest="output",help="output file name")
    p.add_argument('-y','--yy',type=str,dest="yy",help="y axis indicator")
    p.add_argument('-all_geneID','--all_geneID',type=str,dest="all_geneID",help="output file name for all DE gene id")    
    p.add_argument('-overlap_geneID','--overlap_geneID',type=str,dest="overlap_geneID",help="output file name for overlapped gene id")
    p.add_argument('-nonOverlap_geneID','--nonOverlap_geneID',type=str,dest="nonOverlap_geneID",help="output file name for nonOverlapped gene id")
    return p.parse_args()

def getGene(csv,bedChrom,bedHashh):
    '''This function will read the csv file from DE analysis, seperate genes into relocated and not_relcated two categories, write into three files for all genes, relocated genes and not_relocated genes, and report gene number for the three categories'''
    '''This function will return two arrays of FPKM log2-fold changes in the two cell lines for relocated and not_relocated genes'''
    global args
    args=ParseArg()
    all_geneIDout=WriteToFile(args.all_geneID)
    overlap_geneIDout=WriteToFile(args.overlap_geneID)
    nonOverlap_geneIDout=WriteToFile(args.nonOverlap_geneID)
    overlapDomainArray=[]
    nonOverlapDomainArray=[]
    allArray=[]
    for line in ReadFromFile(csv):
        row = line.split(',')
        if row[1].strip('"')=="gene":
            continue        
        chrom="chr"+ str(row[15].strip('"'))
        start=row[16]
        end=row[17]        
        is_overlap = False
        for chr in bedChrom:
            if chr==chrom:
                for region in bedHashh[chrom]:
                    if int(start)>=int(region[0])and int(end)<=int(region[1]): # the entire gene must be within a domain to be called
                        log2FoldChange=row[3].strip('-')
                        geneID=row[13].strip('"')
                        overlapDomainArray.append(float(log2FoldChange))
                        allArray.append(float(log2FoldChange))
                        print >>overlap_geneIDout, '%s\t%s\t%s\t%s' % (str(geneID),chrom,start,end)
                        print >>all_geneIDout, '%s\t%s\t%s\t%s' % (str(geneID),chrom,start,end)
                        is_overlap = True
                        break
        if is_overlap == False:
                log2FoldChange=row[3].strip('-')
                geneID=row[13].strip('"')
                nonOverlapDomainArray.append(float(log2FoldChange))
                allArray.append(float(log2FoldChange))
                print >>nonOverlap_geneIDout, '%s\t%s\t%s\t%s' % (str(geneID),chrom,start,end)
                print >>all_geneIDout, '%s\t%s\t%s\t%s' % (str(geneID),chrom,start,end)
    print "gene (overlap with domain) num:"+str(len(overlapDomainArray))
    print "gene (not overlap with domain) num:"+str(len(nonOverlapDomainArray))
    print "all gene num:"+str(len(allArray))
    return overlapDomainArray,nonOverlapDomainArray

def Main():
    global args
    args=ParseArg()
    '''Read the bed file and collect all the regions by returning a dictionary (chr as keys)'''
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
            array.append((row[1],row[2],row[3]))
        else:
            array.append((row[1],row[2],row[3]))
    bedHash[chrom]=array
    print bedChr
    chrom=None
    '''Get two lists of FPKM log2-fold changes for relocated and not_relocated genes, and calculate P-value by two-tailed t-test'''
    overlap, nonOverlap =getGene(args.csv,bedChr,bedHash) 
    print "overlap vs nonoverlap log2FC #results > (t-statistic, two-tailed p-value)    "+str(stats.ttest_ind(overlap,nonOverlap, equal_var = False))   #results > (t-statistic, two-tailed p-value)
    '''Generate a box plot for the two lists'''
    data=[overlap,nonOverlap]
    sns.set()
    fig,ax=plt.subplots()
    ax=sns.boxplot(data=data)
    plt.xticks([0, 1], ["relocated","not_relocated"])
    plt.ylabel("log2 fold change"+args.yy,fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20, width=1.5)
    plt.subplots_adjust(left=.2)
    plt.savefig(args.output+'.eps',format='eps') 
    plt.close()

if __name__=="__main__":
    Main()
