#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
import seaborn as sns
import pandas as pd
import os
import operator


def main() :
    """Gathers arguments and runs functions"""
    parser = argparse.ArgumentParser(description='Creates a coverage plot.')
    parser.add_argument('cov',nargs=1,type=str,
                        help='File of filenames newline separated coverage files output of COCOA.')
    parser.add_argument('--output','-o',required=False,nargs=1,type=str,default=["Filename"],
                        help='Prefix of output files.')
    parser.add_argument('--interactive',required=False,type=bool,nargs='?',const=True,default=False,
                        help='Interactive plot, zoom, crop, etc.')
    args = parser.parse_args()

    covFoFn = args.cov[0]
    prefix = args.output[0]
    showyes = False
    if args.interactive :
        showyes = True

    FoFn = open(covFoFn, "r")
    covList = {}
    lenList = []
    meanList = []
    names = {}
    i = 0
    for line in FoFn :
        fn = line.strip()
        s = fn.split("\t")
        vCov, length, mean = readCov(s[0])
        covList[i] = vCov
        meanList.append(mean)
        lenList.append(length)
        try :
            names[i] = s[1]
        except :
            names[i] = os.path.basename(s[0])
        i+=1

    keyOrder = []
    cList = []
    nDict = {}
    for k in sorted( covList, key=lambda k: len(covList[k]) ) :
        keyOrder.append(k)
    for i,k in enumerate(keyOrder) :
        cList.append(covList[k])
        nDict[i] = names[k]

    del covList
    del names

    mxLen = max(lenList)
    for Vcov in cList :
        while len(Vcov) != mxLen :
            Vcov.append(np.NaN)

    covDF = pd.DataFrame(cList)
    covDF.rename(index=nDict,inplace=True)
    finalMean = np.mean(meanList)
    doPlot(covDF, prefix, showyes, finalMean)

def readCov(fn) :
    vCov = []
    f = open(fn, "r")
    for line in f :
        s = line.split("\t")
        vCov.append(int(s[2]))
    f.close()
    return vCov, len(vCov), np.mean(vCov)


def doPlot(df, prefix, showyes, mean) :
    print(df)
    plt.style.use("seaborn")
    fig, ax = plt.subplots(figsize=(25,0.5*len(df)))

    c = ax.pcolor(df, cmap="viridis", vmin=0, vmax=2*mean)
    ax.set_yticks( np.arange(df.shape[0])+0.5, minor=False )
    ax.set_yticklabels(df.index.tolist())
    ax.set_xlabel('Position (bp)', fontsize=14)
    fig.colorbar(c, ax=ax)
    if showyes :
        plt.show()
    plt.savefig(prefix + ".png", format="png", dpi=300)

if __name__ == "__main__" :
    main()
