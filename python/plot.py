#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def main() :
    """Gathers arguments and runs functions"""
    parser = argparse.ArgumentParser(description='Creates a coverage plot.')
    parser.add_argument('cov',nargs=1,type=str,
                        help='Coverage file output of COCOA.')
    parser.add_argument('ref',nargs=1,type=str,
                        help='Fasta file containing the contig to plot')
    parser.add_argument('ctg',nargs=1,type=str,
                        help='Contig/Scaffold name (must be in fasta/fastq file)')
    parser.add_argument('--output','-o',required=False,nargs=1,type=str,default=["Filename"],
                        help='Prefix of output files.')
    parser.add_argument('--interactive',required=False,type=bool,nargs='?',const=True,default=False,
                        help='Interactive plot, zoom, crop, etc.')
    args = parser.parse_args()

    cov = args.cov[0]
    ref = args.ref[0]
    ctg = args.ctg[0]
    prefix = args.output[0]
    showyes = False
    if args.interactive :
        showyes = True

    seq = readRef(ref, ctg)
    length = len(seq)
    vCov, start, end = readCov(cov)

    try :
        seqToPlot = seq[start:end+1]
        end += 1
    except : # JUST IN CASE IT FAILS WHEN WHOLE CONTIG IS TO BE PLOTTED
        seqToPlot = seq[start:end]
    del seq
    lengthToPlot = len(seqToPlot)
    del length
    vGC = cGC(seqToPlot, lengthToPlot)
    if len(vCov) != lengthToPlot : # SHOULD NOT HAPPEN
        print("WARNING: length of coverage and length of sequence are different!")
        print("\t", len(vCov), " diff ", lengthToPlot)
    doPlot(lengthToPlot, vGC, vCov, ctg, start, end, prefix, showyes)

def readRef(ref, ctg) :
    seq = ""
    ln = 0
    f = open(ref, "r")
    count = False
    for line in f :
        if ">"+ctg in line :
            count = True
        if count :
            if line[0] != ">" :
                seq += line.strip() # LOADS SEQUENCE IN MEMORY
    f.close()
    return seq

def cGC(seq, length) :
    ratio = []
    wdw = int(0.01*length) # WINDOWS SIZE OF 0.01 length of seq
    for i in range(length-wdw) :
        window = seq[i:i+wdw]
        GC = 0
        total = 0
        for nucl in window :
            if nucl in "GC" :
                GC += 1
            total += 1
        ratio.append( float(GC/total) )
    return ratio

def readCov(cov) :
    start = None
    end = 0
    vCov = []
    f = open(cov, "r")
    for line in f :
        s = line.split("\t")
        if start == None :
            start = int(s[1])
        vCov.append(int(s[2]))
        end = int(s[1])
    f.close()
    return vCov, start, end


def doPlot(lengthToPlot, vGC, vCov, ctg, start, end, prefix, showyes) :

    plt.style.use("seaborn")
    fig, ax1 = plt.subplots(figsize=(15,15))

    ax1.set_title('{} : {} - {}'.format(ctg, start, end), fontsize=20)
    x = [x for x in range(start, end)]
    while (len(vGC) < len(x)) :
        vGC.append(np.NaN)
    ax1.plot(x, vCov, linewidth=0.2, color="red", label="Coverage")
    ax1.set_ylabel('Coverage', fontsize=14)
    ax1.set_xlabel('Position (bp)', fontsize=14)
    ax1.set_ylim(0,)
    ax1.set_xlim(start,end)
    ax1.fill_between(x, vCov, alpha=0.5, color="red")
    ax2 = ax1.twinx()
    ax2.plot(x, vGC, linewidth=0.2, color="black", label="GC")
    ax2.set_ylim(0,1)
    ax2.set_xlim(start,end)
    ax2.set_ylabel("GC")
    ax2.grid(False)
    #ax2.tick_params('y', colors="black")

    if showyes :
        plt.show()
    plt.savefig(prefix + ".png", format="png", dpi=300)

if __name__ == "__main__" :
    main()
