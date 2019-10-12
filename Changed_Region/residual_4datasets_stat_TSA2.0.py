#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 20190306
# Last-modified: 

'''This code will compare 4 wig files (4 replicates from 2 cell lines), calculate residuals with 4 possibilities for each value, fit the residues into gaussian distribution, and return thresholds'''

import os,sys,argparse
from TSA_utility import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from TSA_utility import *
import scipy.stats
from scipy.stats import norm

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-c1r1','--c1r1',type=str,dest="c1r1",help="wig file1 (e.g. TSA-seq percentile wig file1, cell line1 replicate1)")
    p.add_argument('-c1r2','--c1r2',type=str,dest="c1r2",help="wig file1 (e.g. TSA-seq percentile wig file1, cell line1 replicate2)")
    p.add_argument('-c2r1','--c2r1',type=str,dest="c2r1",help="wig file1 (e.g. TSA-seq percentile wig file1, cell line2 replicate1)")
    p.add_argument('-c2r2','--c2r2',type=str,dest="c2r2",help="wig file1 (e.g. TSA-seq percentile wig file1, cell line2 replicate2)")
    p.add_argument('-o','--output',type=str,dest="output",help="output gaussian distribution figure file name")
    p.add_argument('-P','--P-value',type=float,dest="P",help="P value")
    return p.parse_args()

def CollectValue(wigfile):
    '''This function will collect all rescaled TSA-Seq scores (20kb bin)'''
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

def Compare(first_wig,second_wig, third_wig,fourth_wig):
    '''This function will calculate residuals between replicates'''
    global args
    first_hash=CollectValue(first_wig)
    second_hash=CollectValue(second_wig)
    third_hash=CollectValue(third_wig)
    fourth_hash=CollectValue(fourth_wig)
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
    array=[]
    for chrom in chrs:       
        A_list = first_hash[chrom]
        B_list = second_hash[chrom]
        C_list = third_hash[chrom]  
        D_list = fourth_hash[chrom]  
        if not (len(A_list) == len(B_list) and len(A_list) == len(C_list) and len(A_list) == len(D_list)):
            print "The bin number for chrs[i] 2 tsa files are different"
            sys.exit("terminated") 
        for j in range (len(A_list)):
            if A_list[j][0]==B_list[j][0] and A_list[j][0]==C_list[j][0] and A_list[j][0]==D_list[j][0]:
                change1=0.5*((float(A_list[j][1])-float(B_list[j][1]))+(float(C_list[j][1])-float(D_list[j][1])))
                change2=0.5*((float(A_list[j][1])-float(B_list[j][1]))+(float(D_list[j][1])-float(C_list[j][1])))
                change3=0.5*((float(B_list[j][1])-float(A_list[j][1]))+(float(C_list[j][1])-float(D_list[j][1])))
                change4=0.5*((float(B_list[j][1])-float(A_list[j][1]))+(float(D_list[j][1])-float(C_list[j][1])))
                array.append(change1)
                array.append(change2)
                array.append(change3)
                array.append(change4)
    return array

def Main():
    global args
    args=ParseArg()
    Stat=Compare(args.c1r1,args.c1r2, args.c2r1,args.c2r2)   
    print "total number in stat array is:"+str(len(Stat))
    print "max statarray value:"+str(max(Stat))
    print "min statarray value:"+str(min(Stat))
    justStd=np.std(Stat)
    justMean=np.mean(Stat)
    mean,std=norm.fit(Stat)
    print "fitted gaussian mean is:"+str(mean)
    print "fitted gaussian std is:"+str(std)
    print "THE array mean is:"+str(justMean)
    print "THE array std is:"+str(justStd)
    high=scipy.stats.norm(mean,std).ppf(1-args.P)
    low=scipy.stats.norm(mean,std).ppf(args.P)
    print "upper threshold with the input P value is:" + str(high)
    print "lower threshold with the input P value is:"+str(low)
    plt.hist(Stat, bins=np.arange(np.floor(min(Stat))-1, np.ceil(max(Stat))+1, 0.1), facecolor='k')
    plt.xlabel('phi')
    plt.ylabel ('count')
    plt.ylim (0,)
    plt.title('Histogram of stat array')
    plt.savefig(args.output+"_statArray.eps",format="eps")
    plt.close()
    plt.xlabel('phi')
    plt.ylabel ('pdf')
    plt.title('pdf of stat array')
    x=np.linspace(min(Stat),max(Stat),100)
    y=norm.pdf(x,mean,std)
    plt.plot(x,y)
    plt.savefig(args.output+"_statArray_pdf.eps",format="eps")
    plt.close()
    plt.hist(Stat, bins=np.arange(np.floor(min(Stat))-1, np.ceil(max(Stat))+1, 0.1), facecolor='k', normed=True)
    plt.xlabel('phi')
    plt.title('Histogram (normed) of stat array')
    x=np.linspace(min(Stat),max(Stat),100)
    y=norm.pdf(x,mean,std)
    plt.plot(x,y)
    plt.savefig(args.output+"_statArray_normed.eps", format="eps")
    plt.close()
    logging("DONE!!!")

if __name__=="__main__":
    Main()
