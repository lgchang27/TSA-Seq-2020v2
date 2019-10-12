#!/usr/bin/python
# Programmer : Liguo Zhang
# Date: 20190307
# Last-modified: 

'''This code will compare two cell lines with two biological replicates each, identify all 20 kb bins above a input threshold (always cell type 2 - cell type 1).
It will take a second threshold of domain size and return merged segments when the region above the size threshold.
It will output a wig file with all residues, a bed and a merged bed file showing the changed regions.
It will also output a replicate-mean signal bw file for each cell line.
It will also output a bin residue distribution and fit of gaussian distribution.
'''

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
    p.add_argument('-o','--output',type=str,dest="output",help="output tsa compare file name")
    p.add_argument('-c','--change_threshold', type=float, dest='change_threshold',help='set the threshold')
    p.add_argument('-g','--genome',type=str,dest="genome",help="human genome file (two column, first col is chromosome name, second chromosome is chromosome size")
    p.add_argument('-w','--window',type=int,dest="window",help="window size of input wigs")
    p.add_argument('-s','--domain_size',type=int,dest="domain_size",help="threshold for domain size")
    p.add_argument('-cell1','--cell1',type=str,dest="cell1",help="mean replicate signal for cell line1")
    p.add_argument('-cell2','--cell2',type=str,dest="cell2",help="mean replicate signal for cell line2")
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

def Compare(first_wig,second_wig, third_wig,fourth_wig, hgenome, win):
    '''This function will calculate residuals between replicates'''
    global args
    first_hash=CollectValue(first_wig)
    second_hash=CollectValue(second_wig)
    third_hash=CollectValue(third_wig)
    fourth_hash=CollectValue(fourth_wig)
    wigout = WriteToFile(args.output + '.wig')
    bedout = WriteToFile(args.output + '.bed')
    bedout2 = WriteToFile("otherway"+args.output + '.bed')
    wig1=WriteToFile(args.cell1+".wig")
    wig2=WriteToFile(args.cell2+".wig")
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
    compare_hash = {}
    array=[]
    residue_all=[]
    for chrom in chrs:       
        chrom_size = int(LoadGenome(hgenome)[chrom])
        A_list = first_hash[chrom]
        B_list = second_hash[chrom] 
        C_list = third_hash[chrom]
        D_list = fourth_hash[chrom]       
        print >>wigout, "variableStep chrom=%s span=%d" % (chrom,win)
        print >>wig1, "variableStep chrom=%s span=%d" % (chrom,win)
        print >>wig2, "variableStep chrom=%s span=%d" % (chrom,win)
        if not ((len(A_list) == len(B_list) and len(A_list) == len(C_list) and len(A_list) == len(D_list))):
            print "The bin number for chrs[i] 2 tsa files are different"
            sys.exit("terminated")  
        for j in range (len(A_list)):
            if (A_list[j][0]==B_list[j][0] and A_list[j][0]==C_list[j][0] and A_list[j][0]==D_list[j][0]):
                change=0.5*((float(C_list[j][1])+float(D_list[j][1]))-(float(A_list[j][1])+float(B_list[j][1])))
                signal1=0.5*(float(A_list[j][1])+float(B_list[j][1]))
                signal2=0.5*(float(C_list[j][1])+float(D_list[j][1]))
                residue_all.append(change)
                if j== len(A_list) - 1 and (chrom_size-int(A_list[j][0])) < win:
                    print >>wigout, "variableStep chrom=%s span=%d" % (chrom,chrom_size-int(A_list[j][0]))
                    print >>wigout, "%d\t%.6f" % (int(A_list[j][0]),float(change))
                    print >>wig1, "variableStep chrom=%s span=%d" % (chrom,chrom_size-int(A_list[j][0]))
                    print >>wig1, "%d\t%.6f" % (int(A_list[j][0]),float(signal1))
                    print >>wig2, "variableStep chrom=%s span=%d" % (chrom,chrom_size-int(A_list[j][0]))
                    print >>wig2, "%d\t%.6f" % (int(A_list[j][0]),float(signal2))
                else:
                    print >>wigout, "%d\t%.6f" % (int(A_list[j][0]),float(change))
                    print >>wig1, "%d\t%.6f" % (int(A_list[j][0]),float(signal1))
                    print >>wig2, "%d\t%.6f" % (int(A_list[j][0]),float(signal2))          
                if change > args.change_threshold:
                    if j== len(A_list) - 1 and (chrom_size-int(A_list[j][0])) < win:
                        print >>bedout, "%s\t%d\t%d\t%.6f" % (chrom,int(A_list[j][0]),chrom_size,float(change))
                    else:
                        print >>bedout, "%s\t%d\t%d\t%.6f" % (chrom,int(A_list[j][0]),int(A_list[j][0])+win,float(change))
                    array.append((A_list[j][0],change))
                if change < (args.change_threshold*(-1)):
                    if j== len(A_list) - 1 and (chrom_size-int(A_list[j][0])) < win:
                        print >>bedout2, "%s\t%d\t%d\t%.6f" % (chrom,int(A_list[j][0]),chrom_size,float(change))
                    else:
                        print >>bedout2, "%s\t%d\t%d\t%.6f" % (chrom,int(A_list[j][0]),int(A_list[j][0])+win,float(change))
                    array.append((A_list[j][0],change))
        compare_hash[chrom]=array
        array=[]
    wigout.flush()
    wig1.flush()
    wig2.flush()
    bedout.flush()
    bedout2.flush()
    wig2bw = "utilities/wigToBigWig -clip %s %s %s" % ( args.output+'.wig', args.genome, args.output+'.bw')
    os.system(wig2bw) 
    wig2bw = "utilities/wigToBigWig -clip %s %s %s" % ( args.cell1+'.wig', args.genome, args.cell1+'.bw')
    os.system(wig2bw) 
    wig2bw = "utilities/wigToBigWig -clip %s %s %s" % ( args.cell2+'.wig', args.genome, args.cell2+'.bw')
    os.system(wig2bw)   
    return compare_hash,residue_all

def MergeSegment_adjacent (filename,outname):
    '''This funciton will merge adjacent bins to call a region above the input threshold'''
    global args
    out = WriteToFile(outname)
    chrom = None
    sizeList=[]
    rsdList=[]
    for line in ReadFromFile(filename):
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        if row[0] != chrom:
            if chrom is not None:
              if(int(stop)-int(start))>=args.domain_size:
                print >>out, "%s\t%s\t%s\t%.6f" % (chrom,start,stop,np.mean(residualList))
                sizeList.append(int(stop)-int(start))
                rsdList.append(abs(np.mean(residualList)))
            chrom = row[0]
            start = row[1]
            stop = row[2]
            residualList=[float(row[3])]
        else:
            if row[1] != stop:
                if(int(stop)-int(start))>=args.domain_size:
                    print >>out, "%s\t%s\t%s\t%.6f" % (chrom,start,stop,np.mean(residualList))
                    sizeList.append(int(stop)-int(start))
                    rsdList.append(abs(np.mean(residualList)))
                chrom = row[0]
                start = row[1]
                stop = row[2]
                residualList=[float(row[3])]
            elif row[1] == stop:
                stop = row[2]
                residualList.append(float(row[3]))
    if (int(stop)-int(start))>=args.domain_size:
        print >>out, "%s\t%s\t%s\t%.6f" % (chrom,start,stop,np.mean(residualList))
        sizeList.append(int(stop)-int(start))
        rsdList.append(abs(np.mean(residualList)))
    out.flush()
    return sizeList,rsdList

def Main():
    global args
    args=ParseArg()
    genome = LoadGenome(args.genome)
    TSAcompare,rsd=Compare(args.c1r1,args.c1r2, args.c2r1,args.c2r2, args.genome,args.window)
    chrs = TSAcompare.keys()
    print 'changed region: chr num:'+str(len(chrs))
    
    '''Identify bins with residual above threshold'''
    TSA_array=[]
    TSA_abs_array=[]
    for i in range(len(chrs)):
        TSA_list = TSAcompare[chrs[i]]      
        for j in range (len(TSA_list)):
                TSA_array.append(TSA_list[j][1])
                TSA_abs_array.append(abs(TSA_list[j][1]))
    print "changed TSA bin num:"+str(len(TSA_array))
    print "max TSA residual above threshold:"+str(max(TSA_array))
    print "min TSA residual above threshold:"+str(min(TSA_array))
    print "min TSA abs residual above threshold:"+str(min(TSA_abs_array))
    plt.hist(TSA_array, bins=np.arange(np.floor(min(TSA_array))-1, np.ceil(max(TSA_array))+1, 1), facecolor='k')
    plt.xlabel('TSA_residual_above_threshold (bin)')
    plt.ylabel ('count')
    plt.ylim (0,)
    plt.title('Histogram of TSA_residule_above_threshold (bin)')
    plt.savefig(args.output+"_bin_hst.png")
    plt.close()

    '''All bin residuals'''
    print "total number of all rsd is:"+str(len(rsd))
    print "max rsd value:"+str(max(rsd))
    print "min rsd value:"+str(min(rsd))
    mean,std=norm.fit(rsd)
    print "fitted gaussian mean is:"+str(mean)
    print "fitted gaussian std is:"+str(std)
    plt.hist(rsd, bins=np.arange(np.floor(min(rsd))-1, np.ceil(max(rsd))+1, 1), facecolor='k')
    plt.xlabel('residue')
    plt.ylabel ('count')
    plt.ylim (0,)
    plt.title('Histogram of residues')
    plt.savefig(args.output+"_AllBinResidues.png")
    plt.close()
    plt.xlabel('residue')
    plt.ylabel ('pdf')
    plt.title('pdf of residues')
    x=np.linspace(min(rsd),max(rsd),100)
    y=norm.pdf(x,mean,std)
    plt.plot(x,y)
    plt.savefig(args.output+"_AllBinResidues_pdf.png")
    plt.close()


    '''One way compare: cell type 2 signal > cell type 1 signal'''
    SortBed(args.output + '.bed')
    mergename = args.output + "_mergeAdjacent.bed"
    size1,rsd1= MergeSegment_adjacent(args.output+'.bed',mergename)
    print "one way size array length:"+str(len(size1))
    print "one way rsd array length:"+str(len(rsd1))

    '''The other way compare: cell type 1 signal > cell type 2 signal'''
    SortBed('otherway'+args.output + '.bed')
    mergename = 'otherway'+args.output + "_mergeAdjacent.bed"
    size2,rsd2= MergeSegment_adjacent('otherway'+args.output+'.bed',mergename)
    print "the other way size array length:"+str(len(size2))
    print "the other way rsd array length:"+str(len(rsd2))

    '''print domain info'''
    print 'one way domain number:'+str(len(size1))+'\nmean domian size:'+str(int(np.mean(size1)))+'\nmedian domian size:'+str(int(np.median(size1)))+'\nmax domian size:'+str(int(np.max(size1)))+'\nmin domian size:'+str(int(np.min(size1)))
    print 'one way domain mean rsd:'+str(float(np.mean(rsd1)))+'\nmedian rsd:'+str(float(np.median(rsd1)))+'\nmax rsd:'+str(float(np.max(rsd1)))+'\nmin rsd:'+str(float(np.min(rsd1)))
    print 'the other way domain number:'+str(len(size2))+'\nmean domian size:'+str(int(np.mean(size2)))+'\nmedian domian size:'+str(int(np.median(size2)))+'\nmax domian size:'+str(int(np.max(size2)))+'\nmin domian size:'+str(int(np.min(size2)))
    print 'the othe way domain mean rsd:'+str(float(np.mean(rsd2)))+'\nmedian rsd:'+str(float(np.median(rsd2)))+'\nmax rsd:'+str(float(np.max(rsd2)))+'\nmin rsd:'+str(float(np.min(rsd2)))

    '''Draw histogram for one way compare: cell type 2 signal > cell type 1 signal'''
    plt.hist(size1, bins=np.arange(np.floor(min(size1))-20000, np.ceil(max(size1))+20000, 20000), facecolor='k')
    plt.xlabel('TSA_residual_above_threshold_domain_size')
    plt.ylabel ('count')
    plt.title('Histogram of TSA_residule_above_threshold_domain_size')
    plt.annotate('domain number:'+str(len(size1))+'\nmean domian size:'+str(int(np.mean(size1)))+'\nmedian domian size:'+str(int(np.median(size1)))+'\nmax domian size:'+str(int(np.max(size1)))+'\nmin domian size:'+str(int(np.min(size1))),xy=(0.5, 0.6), xycoords='axes fraction',fontsize=10)
    plt.savefig(args.output+"_domain_size_hst.png")
    plt.close()
    sort="sort -k1,1 -k2,2n %s > %s" % (args.output+"_mergeAdjacent.bed",args.output+"_mergeAdjacent.sorted.bed")
    bed2bb = "utilities/bedToBigBed %s %s %s" % ( args.output+'_mergeAdjacent.sorted.bed', args.genome, args.output+'_mergeAdjacent.bb')
    os.system(sort)
    os.system(bed2bb)
    plt.hist(rsd1, bins=np.arange(np.floor(min(rsd1))-1, np.ceil(max(rsd1))+1, 1), facecolor='k')
    plt.xlabel('TSA_residual_above_threshold_domain')
    plt.ylabel ('count')
    plt.title('Histogram of TSA_residule_above_threshold_domain')
    plt.annotate('mean rsd:'+str(float(np.mean(rsd1)))+'\nmedian rsd:'+str(float(np.median(rsd1)))+'\nmax rsd:'+str(float(np.max(rsd1)))+'\nmin rsd:'+str(float(np.min(rsd1))),xy=(0.5, 0.6), xycoords='axes fraction',fontsize=10)
    plt.savefig(args.output+"_domain_hst.png")
    plt.close()    

    '''Draw histogram for the other way compare: cell type 1 signal > cell type 2 signal'''
    plt.hist(size2, bins=np.arange(np.floor(min(size2))-20000, np.ceil(max(size2))+20000, 20000), facecolor='k')
    plt.xlabel('TSA_residual_above_threshold_domain_size')
    plt.ylabel ('count')
    plt.title('Histogram of TSA_residule_above_threshold_domain_size')
    plt.annotate('domain number:'+str(len(size2))+'\nmean domian size:'+str(int(np.mean(size2)))+'\nmedian domian size:'+str(int(np.median(size2)))+'\nmax domian size:'+str(int(np.max(size2)))+'\nmin domian size:'+str(int(np.min(size2))),xy=(0.5, 0.6), xycoords='axes fraction',fontsize=10)
    plt.savefig('otherway'+args.output+"_domain_size_hst.png")
    plt.close()
    sort2="sort -k1,1 -k2,2n %s > %s" % ('otherway'+args.output+"_mergeAdjacent.bed",'otherway'+args.output+"_mergeAdjacent.sorted.bed")
    bed2bb2 = "utilities/bedToBigBed %s %s %s" % ( 'otherway'+args.output+'_mergeAdjacent.sorted.bed', args.genome, 'otherway'+args.output+'_mergeAdjacent.bb')
    os.system(sort2)
    os.system(bed2bb2)
    plt.hist(rsd2, bins=np.arange(np.floor(min(rsd2))-1, np.ceil(max(rsd2))+1, 1), facecolor='k')
    plt.xlabel('TSA_residual_above_threshold_domain')
    plt.ylabel ('count')
    plt.title('Histogram of TSA_residule_above_threshold_domain')
    plt.annotate('mean rsd:'+str(float(np.mean(rsd2)))+'\nmedian rsd:'+str(float(np.median(rsd2)))+'\nmax rsd:'+str(float(np.max(rsd2)))+'\nmin rsd:'+str(float(np.min(rsd2))),xy=(0.5, 0.6), xycoords='axes fraction',fontsize=10)
    plt.savefig('otherway'+args.output+"_domain_hst.png")
    plt.close()  

    logging("DONE!!!")

if __name__=="__main__":
    Main()
