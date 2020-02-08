#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 20191215

import os,sys,argparse
from TSA_utility import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-b1','--bed1',type=str,dest="bed1",help="bed file 1")
    p.add_argument('-b2','--bed2',type=str,dest="bed2",help="bed file 2")
    p.add_argument('-o','--output',type=str,dest="output",help="output file name")
    return p.parse_args()

def Main():
    global args
    args=ParseArg()
    array1=[]
    array2=[]
    for line in ReadFromFile(args.bed1):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        array1.append(abs(float(row[5])))
    print array1
    for line in ReadFromFile(args.bed2):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        array2.append(abs(float(row[5])))
    print array2

    array=array1+array2
    n=len(array)
    maxd = max(array)
    mind = min(array)
    mean = np.mean(array)
    median = np.median(array)
    std = np.std(array)
    print "n:"+str(n)
    print "max:"+str(maxd)
    print "min:"+str(mind)
    print "mean:"+str(mean)
    print "median:"+str(median)
    print "std:"+str(std)

    plt.hist(array, bins=np.arange(0.1, 0.35, 0.01), facecolor='k')
    plt.xlabel('distance residual (um)', fontsize=20)
    plt.ylabel ('count',fontsize=20)
    plt.annotate('max:'+str(np.around(maxd,decimals=3))+'\nmin:'+str(np.around(mind,decimals=3))+'\nmean:'+str(np.around(mean,decimals=3))+'\nmedian:'+str(np.around(median,decimals=3)),xy=(0.5, 0.6), xycoords='axes fraction',fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20, width=2)
    plt.tick_params(axis='both', which='minor', labelsize=20, width=2)
    #plt.subplot().set_xticks([0,0.1,0.2,0.3])
    plt.subplot().spines['right'].set_visible(False)
    plt.subplot().spines['top'].set_visible(False)
    plt.subplot().spines['left'].set_linewidth(2)
    plt.subplot().spines['bottom'].set_linewidth(2)    
    plt.subplots_adjust(bottom=.2, left=.2)
    plt.savefig(args.output+"_hist.eps", format='eps')
    plt.close()

    logging("DONE!!!")

if __name__=="__main__":
    Main()
