#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 
# Last-modified: 27 Sep 2019   

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
    p.add_argument('-w1','--wig1',type=str,dest="wig1",help="percentile wig file1 (e.g. wig file1)")
    p.add_argument('-w2','--wig2',type=str,dest="wig2",help="percentile wig file2 (e.g. wig file2)")    
    p.add_argument('-o','--output',type=str,dest="output",help="output file name")
    p.add_argument('-x','--x',type=str,dest="x",help="x axis name")
    p.add_argument('-y','--y',type=str,dest="y",help="y axis name")
    return p.parse_args()

def CollectValue(wigfile):
    '''This funciton collected all percentile values (20kb bin) from wig file'''
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
            logging("Process: %d line processed in wig1" % (num))
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
    TSApercentile1=CollectValue(args.wig1)
    TSApercentile2=CollectValue(args.wig2)
    chrs_1 = TSApercentile1.keys()
    chrs_2 = TSApercentile2.keys()
    chrs = list(set(chrs_1) & set(chrs_2))
    print chrs
    print 'chr num:'+str(len(chrs))
    TSApercentile1_array=[]
    TSApercentile2_array=[]
    for i in range(len(chrs)):
        A_list = TSApercentile1[chrs[i]]
        B_list = TSApercentile2[chrs[i]]       
        if len(A_list) != len(B_list):
            print "The bin number for chrs[i] in the two tsa files are different"
            sys.exit("terminated") 
        for j in range (len(A_list)):
            if A_list[j][0]==B_list[j][0]:
                TSApercentile1_array.append(float(A_list[j][1])+50.0)
                TSApercentile2_array.append(float(B_list[j][1])+50.0)
    print max(TSApercentile1_array)
    print min(TSApercentile1_array)
    print max(TSApercentile2_array)
    print min(TSApercentile2_array)
    print "TSApercentile1 bin num:"+str(len(TSApercentile1_array))
    print "TSApercentile2 bin num:"+str(len(TSApercentile2_array))
    
    plt.figure(figsize=(7, 7))
    x=TSApercentile1_array
    y=TSApercentile2_array
    plt.hist2d(x,y,(200,200),norm=matplotlib.colors.LogNorm(),cmap=plt.cm.jet)
    plt.colorbar()
    plt.axes().set_aspect('equal')
    plt.xlabel(args.x,fontsize=30)
    plt.ylabel(args.y,fontsize=30)
    plt.title('SON TSA Percentile',size=30)
    plt.xlim (-5,105)
    plt.ylim (-5,105)
    plt.tick_params(axis='both', which='major', labelsize=30,width=2)
    plt.axes().set_aspect('equal')
    plt.subplot().set_xticks([0,20,40,60,80,100])
    plt.subplot().set_yticks([0,20,40,60,80,100])
    x1,y1=[90,90],[100,90]
    x2,y2=[100,90],[90,90]
    x3,y3=[90,100],[100,100]
    x4,y4=[100,100],[100,90]
    plt.subplot().spines['left'].set_linewidth(2)
    plt.subplot().spines['bottom'].set_linewidth(2)
    plt.subplot().spines['right'].set_linewidth(2)
    plt.subplot().spines['top'].set_linewidth(2)
    plt.subplot().spines['right'].set_visible(False)
    plt.subplot().spines['top'].set_visible(False)
    plt.subplots_adjust(bottom=.18, left=.2)
    plt.plot(x1,y1,linewidth=3.5,linestyle='--',color='red')
    plt.plot(x2,y2,linewidth=3.5,linestyle='--',color='red') 
    plt.plot(x3,y3,linewidth=3.5,linestyle='--',color='red')
    plt.plot(x4,y4,linewidth=3.5,linestyle='--',color='red')
    plt.savefig(args.output,format='eps')
    plt.close()
    
if __name__=="__main__":
    Main()
