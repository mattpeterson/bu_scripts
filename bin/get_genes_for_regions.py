#!/usr/bin/env python
import sys

from transnet.chipseq import poisson
import transnet.genome as genome

if __name__ == "__main__":
    file_handle = open(sys.argv[1])
    annotation_handle = open(sys.argv[2])
    annotation = genome.read(annotation_handle, "broad")

    for peak in poisson.parse(file_handle, "1"):
        features = peak.get_regulated_genes(annotation)
        print " ".join([str(f) for f in features])
