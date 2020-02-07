#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 
# Last-modified: 01 Mar 2019

import os,sys,argparse
import numpy as np
import tempfile
from TSA_utility import *
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-w1','--wig1',type=str,dest="wig1",help="TSA-Seq score wig file1 (e.g. wig file1)")
    p.add_argument('-w2','--wig2',type=str,dest="wig2",help="TSA-Seq score wig file2 (e.g. wig file2)")    
    p.add_argument('-x','--x',type=str,dest="xaxis",help="x name")
    p.add_argument('-y','--y',type=str,dest="yaxis",help="y name")
    p.add_argument('-o','--out',type=str,dest="out",help="output name")
    return p.parse_args()

def CollectValue(wigfile):
    '''This funciton collected all TSA-Seq scores (20kb bin) from wig file'''
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
    TSAscore1=CollectValue(args.wig1)
    TSAscore2=CollectValue(args.wig2)
    chrs_1 = TSAscore1.keys()
    chrs_2 = TSAscore2.keys()
    chrs = list(set(chrs_1) & set(chrs_2))
    print chrs
    print 'chr num:'+str(len(chrs))
    TSAscore1_array=[]
    TSAscore2_array=[]
    for i in range(len(chrs)):
        A_list = TSAscore1[chrs[i]]
        B_list = TSAscore2[chrs[i]]       
        if len(A_list) != len(B_list):
            print "The bin number for chrs[i] in the two tsa files are different"
            sys.exit("terminated") 
        for j in range (len(A_list)):
            if A_list[j][0]==B_list[j][0]:
                TSAscore1_array.append(float(A_list[j][1]))
                TSAscore2_array.append(float(B_list[j][1]))
    print max(TSAscore1_array)
    print min(TSAscore1_array)
    print max(TSAscore2_array)
    print min(TSAscore2_array)
    print "TSAscore1 bin num:"+str(len(TSAscore1_array))
    print "TSAscore2 bin num:"+str(len(TSAscore2_array))
    
    pearR=np.corrcoef(TSAscore1_array,TSAscore2_array)[1,0]
    print pearR

    plt.figure(figsize=(7, 7))
    x=TSAscore1_array
    y=TSAscore2_array
    plt.hist2d(x,y,(200,200),[[-2, 4], [-2.5, 4]],norm=matplotlib.colors.LogNorm(),cmap=plt.cm.jet)
    plt.colorbar()
    plt.axes().set_aspect('equal')
    plt.xlabel(args.xaxis,fontsize=30)
    plt.ylabel(args.yaxis,fontsize=30)
    plt.tick_params(axis='both', which='major', labelsize=30,width=2)
    plt.subplot().set_yticks([-2,-1,0,1,2,3,4])
    plt.subplot().set_xticks([-2,-1,0,1,2,3,4])
    plt.subplot().spines['left'].set_linewidth(2)
    plt.subplot().spines['bottom'].set_linewidth(2)
    plt.subplot().spines['right'].set_linewidth(2)
    plt.subplot().spines['top'].set_linewidth(2)
    plt.subplot().spines['right'].set_visible(False)
    plt.subplot().spines['top'].set_visible(False)
    plt.subplots_adjust(bottom=.2, left=.2)
    plt.savefig(args.out+'.eps',format='eps')
    plt.close()

    
if __name__=="__main__":
    Main()
