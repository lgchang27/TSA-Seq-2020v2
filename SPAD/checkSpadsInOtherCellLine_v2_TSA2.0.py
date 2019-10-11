#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 20180424
# Last-modified: 20181120, optimized wig file reading method
# Check percentile distribution of domains (SPADs) in other cell lines.

import os,sys,argparse
from TSA_utility import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-b','--bed',type=str,dest="bed",help="bed file (e.g. merged bed file for SPADs identification)")    
    p.add_argument('-w1','--wig1',type=str,dest="wig1",help="wig file (e.g. TSA-seq percentile wig file)")
    p.add_argument('-w2','--wig2',type=str,dest="wig2",help="wig file (e.g. TSA-seq percentile wig file)")
    p.add_argument('-w3','--wig3',type=str,dest="wig3",help="wig file (e.g. TSA-seq percentile wig file)")
    p.add_argument('-o','--output',type=str,dest="output",help="output tsa compare file name")
    p.add_argument('-g','--genome',type=str,dest="genome",help="human genome file (two column, first col is chromosome name, second chromosome is chromosome size")
    p.add_argument('-win','--window',type=int,dest="window",help="window size of input wigs")
    return p.parse_args()

def getDomain(bedf):
    '''This function will read a bed file and return all regions'''
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
            array.append((row[1],row[2],row[3]))
        else:
            array.append((row[1],row[2],row[3]))
    bedHash[chrom]=array
    print "bed chr:"
    print bedChr
    print bedHash
    chrom=None
    return bedChr, bedHash

def getDomainMeanPercentiles(chrs,regions,wigf):
    '''This function will collect all regions from the input and calculate region mean percentile in a different cell line with a wig file'''
    wighash={}
    array_temp=[]
    for line in ReadFromFile(wigf):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        if row[0] == 'track':
            continue
        if row[0] == 'variableStep': # new track chrom begin
            if len(array_temp)!=0:
                if chrom != row[1].split('=')[1]:
                    wighash[chrom]=array_temp
                    array_temp=[]
            chrom = row[1].split('=')[1]
        else:
            array_temp.append((float(row[0]),float(row[1])))
            wighash[chrom]=array_temp
    domainMeanPercentileList=[]
    for i in range (len(chrs)):
        domainChr=chrs[i]
        domainList=regions[domainChr]
        for j in range (len(domainList)):
            domainStart=float(domainList[j][0])
            domainEnd=float(domainList[j][1])
            tempList=[]
            if domainChr in wighash.keys():  
                chrList=wighash[domainChr]
                for mm in range (len(chrList)):
                    wigBp=float(chrList[mm][0])
                    wigPercentile=float(chrList[mm][1])+50
                    if wigBp>=domainStart and wigBp<domainEnd:
                            tempList.append(wigPercentile)
                meanPercentile=np.mean(tempList)
                domainMeanPercentileList.append(meanPercentile)
    return domainMeanPercentileList

def Main():
    global args
    args=ParseArg()
    domainChr,domainHash=getDomain(args.bed)
    list1=getDomainMeanPercentiles(domainChr,domainHash,args.wig1)   
    list2=getDomainMeanPercentiles(domainChr,domainHash,args.wig2)   
    list3=getDomainMeanPercentiles(domainChr,domainHash,args.wig3)
    print "domain number check1:"+str(len(list1))
    print "domain number check2:"+str(len(list2))
    print "domain number check3:"+str(len(list3))
    data=[list1,list2,list3]
    sns.set()
    fig,ax=plt.subplots()
    ax=sns.boxplot(data=data)
    ax=sns.swarmplot(data=data,color="0.4")
    plt.xticks([0, 1, 2], ["cell1", "cell2","cell3"]) # change accordingly
    plt.ylim(80,100)
    plt.ylabel('Region mean percentile',fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=20, width=1.5)
    plt.subplots_adjust(left=.2)
    plt.savefig(args.output,format='eps') 
    plt.close()
    logging("DONE!!!")

if __name__=="__main__":
    Main()
