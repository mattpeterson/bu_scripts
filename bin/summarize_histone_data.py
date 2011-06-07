#!/usr/bin/env python

import sys
from glob import glob
import os.path

def read_rpkm_vals(file_to_read):

    symbols = {}
    rpkms = {}

    with open(file_to_read) as file_handle:
        for line in file_handle:
            tokens = line.strip("\r\n").split("\t")
            if len(tokens) == 1:
                continue
            symbols[tokens[0]] = tokens[1]
            rpkms[tokens[0]] = float(tokens[2])

    return symbols, rpkms

def read_scores(scores_file):
    scores = {}
    with open(scores_file) as file_handle:
        for line in file_handle:
            tokens = line.strip("\r\n").split("\t")
            scores[tokens[0]] = float(tokens[1])
    return scores

if __name__ == "__main__":

    files = glob(sys.argv[1] + "/*.scores")
    symbols, rnap = read_rpkm_vals(sys.argv[2])
    symbols, rnaseq = read_rpkm_vals(sys.argv[3])

    values = {}

    for f in files:
        name, ext = os.path.splitext(os.path.basename(f))
        values[name] = read_scores(f)

    print "Gene\tSymbol\tRNAP (RPKM)\tRNASeq(RPKM)\t%s" % ("\t".join(values.keys()))
    for k in rnap.keys():
        out_vals = []
        symbol = symbols[k]
        for k2 in values.keys():
            if k in values[k2]:
                out_vals.append(str(values[k2][k]))
            else:
                out_vals.append(str("0.0"))
        print "%s\t%s\t%f\t%f\t%s" % (k, symbol, rnap[k], rnaseq[k], "\t".join(out_vals))

