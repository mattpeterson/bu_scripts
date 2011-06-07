#!/usr/bin/env python

import sys
from transnet import chipseq, genome

import neurospora_utils

if __name__ == "__main__":
    histone_handle = open(sys.argv[1])
    annotation_handle = open(sys.argv[2])

    annotation = genome.read(annotation_handle, format = "broad", mapping =
                             neurospora_utils.get_contig_mapping())

    peaks = chipseq.sicer.parse(histone_handle, "sicer_rb")
    scores = {}
    for p in peaks:
        features = p.get_regulated_genes(annotation)
        for f in features:
            if f.__class__ == genome.IntergenicRegion:
                if f.left_gene in scores:
                    if p.score > scores[f.left_gene]:
                        scores[f.left_gene] = p.score
                else:
                    scores[f.left_gene] = p.score

                if f.right_gene in scores:
                    if p.score > scores[f.right_gene]:
                        scores[f.right_gene] = p.score
                else:
                    scores[f.right_gene] = p.score
            else:
                if f in scores:
                    if p.score > scores[f]:
                        scores[f] = p.score
                else:
                    scores[f] = p.score

    for k,v in scores.items():
        print("%s\t%f" % (k.locus, v))

