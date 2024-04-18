#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re, os
import pandas as pd
from multiprocessing import Pool

def getOneColumnToList(targetFile, sep="\t", col=0):
    df = pd.read_csv(targetFile, sep=sep)
    return df.iloc[:, col].to_list()

def validateFile(myFile):
    if not os.path.exists(myFile):
        raise Exception("File '%s' not found! Please input again!" % myFile)

    if not os.path.isfile(myFile):
        raise Exception("File '%s' is not a file! Please input again!" % myFile)

    return True

def validateDir(myDir):
    if not os.path.exists(myDir):
        raise Exception("Dir '%s' not found! Please input again!" % myDir)

    if not os.path.isdir(myDir):
        raise Exception("Dir '%s' is not a directory! Please input again!" % myDir)

    return True

def resolveDir(dirName, chdir=True):
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    if chdir:
        os.chdir(dirName)


def parallel_plot(myFunc, values, threads=10):
    # for i in values:
    #     print(i)
    # try:
    #     p = Pool(processes=threads)
    #     results = [p.apply_async(myFunc, param) for param in values]
    #     p.close()
    #     p.join()
    # except KeyboardInterrupt as e:
    #     p.terminate()
    #     p.join()
    #     raise e
    # [print(*param, len(param)) for param in values]
    p = Pool(processes=threads)
    results = [p.apply_async(myFunc, param) for param in values]
    p.close()
    p.join()

    return [res.get() for res in results]

    # for i in values:
    #     myFunc(*i)

def ped2fasta(ped_file, out_file):
    out = open(out_file, 'w')
    with open(ped_file) as f:
        for line in f.readlines():
            tmp = line.strip().split(' ')
            name = tmp[1]
            seq = ''.join(tmp[6:])
            seq = re.sub("0", "N", seq)
            out.write('>'+name+"\n"+seq+"\n")
    out.close()