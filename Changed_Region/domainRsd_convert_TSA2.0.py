#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 20190605

'''
This code will take the _mergeAdjacent.bed file for the changed domains, get the region mean percentile/maxminNormScore in both cell lines (replicate mean) for each region, and convert the percentile/maxminNormScore to distance based on a provided score_distance dictionary in K562
This code will calculate the distance residuals between the converted distences in the two cell lines for each region, always -(cell line 2 - cell line 1).
This code will output a .bed file with five columns: chrom, start, end, distance(cell1), distance(cell2), distance residual (cell1-cell2)
'''

import os,sys,argparse
from TSA_utility import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-b','--bed',type=str,dest="bed",help="bed file ")
    p.add_argument('-d','--dic',type=str,dest="dic",help="dic file ")
    p.add_argument('-p1','--perc1',type=str,dest="perc1",help="TSA-Seq percentile or miaxmin normalized score file1")
    p.add_argument('-p2','--perc2',type=str,dest="perc2",help="TSA-Seq percentile or miaxmin normalized score file2")
    p.add_argument('-p3','--perc3',type=str,dest="perc3",help="TSA-Seq percentile or miaxmin normalized score file3")
    p.add_argument('-p4','--perc4',type=str,dest="perc4",help="TSA-Seq percentile or miaxmin normalized score file4")
    p.add_argument('-o','--output',type=str,dest="output",help="output file name")
    return p.parse_args()

def CollectValue(wigfile):
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
    print "all value num:"+str(len(all))
    print "max vallue:"+str(max(all))
    print "min vallue:"+str(min(all))
    return hash


def Main():
    global args
    args=ParseArg()
    perc_hash1 = CollectValue (args.perc1)
    perc_hash2 = CollectValue (args.perc2)
    perc_hash3 = CollectValue (args.perc3)
    perc_hash4 = CollectValue (args.perc4)
    
    '''read dic'''
    dic={}
    for line in ReadFromFile(args.dic):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        if row[0] == 'percentile':
            continue
        dic[str(row[0])]=float(row[2])   # maxmin score and predicted mean distance in K562
    print dic

    newBed=WriteToFile(args.output)
    '''read bed file for the changed regions'''
    chrom=None
    for line in ReadFromFile(args.bed):
        row = line.strip().split()
        chrom = row[0]
        start = row[1]
        end = row[2]
        perc_list1 = perc_hash1 [chrom]
        perc_list2 = perc_hash2 [chrom]
        perc_list3 = perc_hash3 [chrom]
        perc_list4 = perc_hash4 [chrom]
        array1 = []
        array2 = []
        array3 = []
        array4 = []
        for i in range (len(perc_list1)):
            if int(perc_list1[i][0]) >= int(start) and int(perc_list1[i][0]) < int(end):
                array1.append(float(perc_list1[i][1]))
        #print array1
        for j in range (len(perc_list2)):
            if int(perc_list2[j][0]) >= int(start) and int(perc_list2[j][0]) < int(end):
                array2.append(float(perc_list2[j][1]))
        #print array2
        for m in range (len(perc_list3)):
            if int(perc_list3[m][0]) >= int(start) and int(perc_list3[m][0]) < int(end):
                array3.append(float(perc_list3[m][1]))
        #print array3
        for n in range (len(perc_list4)):
            if int(perc_list4[n][0]) >= int(start) and int(perc_list4[n][0]) < int(end):
                array4.append(float(perc_list4[n][1]))
        #print array4
        if len(array1) != len(array2):
            print "error"
        if len(array1) != len(array3):
            print "error"
        if len(array1) != len(array4):
            print "error"
        cell1rep1 = np.mean(array1)+50
        cell1rep2 = np.mean(array2)+50
        cell2rep1 = np.mean(array3)+50
        cell2rep2 = np.mean(array4)+50
        cell1 = (cell1rep1+cell1rep2)/2
        cell2 = (cell2rep1+cell2rep2)/2
        
        distance1=dic[str(int(np.round(cell1)))]
        distance2=dic[str(int(np.round(cell2)))]   
        distanceRsd = -(distance2 - distance1)        ######for distance rsd, take the negate
        print >>newBed, "%s\t%s\t%s\t%f\t%f\t%f" % (chrom, start, end, distance1, distance2, distanceRsd)          

if __name__=="__main__":
    Main()
