#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 20190608

'''This code will take a gene list, according to coordinates to find the mean TSA-Seq values in two cell lines, and draw a dot plot, x is cell1 and y is cell2'''
'''Can be used to check TSA-Seq scores, percentile, or maxmin rescaled scores'''

import os,sys,argparse
from TSA_utility import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from bx.bbi.bigwig_file import BigWigFile


def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-p1','--p1',type=str,dest="percentile1",help="percentile file1")
    p.add_argument('-p2','--p2',type=str,dest="percentile2",help="percentile file2")
    p.add_argument('-g','--gene',type=str,dest="geneList",help="gene list file")
    p.add_argument('-o','--output',type=str,dest="output",help="output file name without .xx")
    p.add_argument('-y','--yy',type=str,dest="y",help="y axis indicator")
    p.add_argument('-x','--xx',type=str,dest="x",help="x axis indicator")
    return p.parse_args()


def Main():
    global args
    args=ParseArg()
    bw1 = BigWigFile(open(args.percentile1))
    bw2 = BigWigFile(open(args.percentile2))
    gout = WriteToFile(args.output + ".list")
    perc_array1=[]
    perc_array2=[]
    for line in ReadFromFile(args.geneList):
        row = line.strip().split()
        gene = row[0]
        chrom = row[1]
        start = int(row[2])
        end = int(row[3])
        array1 = bw1.get_as_array(chrom,start,end)
        array2 = bw2.get_as_array(chrom,start,end)
        if array1 is not None and array2 is not None:
            perc1 = np.mean(array1)+50
            perc2 = np.mean(array2)+50
            print >>gout, '%s\t%s\t%d\t%d\t%f\t%f' % (gene,chrom,start,end,perc1,perc2)
            perc_array1.append(perc1)
            perc_array2.append(perc2)

    '''scatter plot'''
    sns.set()
    plt.scatter(perc_array1, perc_array2, marker=',', color='black', s=1,alpha=0.1)
    plt.axes().set_aspect('equal')
    plt.xlabel(args.x,fontsize=20)
    plt.ylabel(args.y,fontsize=20)
    plt.ylim (0,100)
    plt.xlim (0,100)
    plt.tick_params(axis='both', which='major', labelsize=20,width=2)
    plt.gca().set_yticks([0,20,40,60,80,100])
    plt.gca().set_xticks([0,20,40,60,80,100])
    x1,y1=[0,89.5],[10.5,100]  
    x2,y2=[10.5,100],[0,89.5] # draw lines showing the threshold to call changed domains
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)
    plt.gca().spines['top'].set_linewidth(2)
    plt.subplots_adjust(bottom=.2, left=.2)
    plt.plot(x1,y1,linewidth=0.5,linestyle='--',color='red')
    plt.plot(x2,y2,linewidth=0.5,linestyle='--',color='red') 
    plt.savefig(args.output+'_dot.eps',format='eps')
    plt.close()

    logging("DONE!!!")

if __name__=="__main__":
    Main()
