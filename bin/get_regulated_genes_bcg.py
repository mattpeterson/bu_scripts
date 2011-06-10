#!/usr/bin/env python

import sys
from transnet import genome
from transnet.chipseq import poisson

def split_features(regulated_genes):
    """
    Go through and classify all of the genes in the region
    """
    genic = []
    upstream = []
    downstream = []

    for r in regulated_genes:
        if r.__class__.__name__ == "Gene":
            genic.append(r.locus)
        elif r.__class__.__name__ == "IntergenicRegion":
            if r.left_gene.strand == "-":
                upstream.append(r.left_gene.locus)
            else:
                downstream.append(r.left_gene.locus)

            if r.right_gene.strand == "+":
                upstream.append(r.right_gene.locus)
            else:
                downstream.append(r.right_gene.locus)

        else:
            # We should never get here
            assert False

    return genic, upstream, downstream

if __name__ == "__main__":
    peaks_handle = open(sys.argv[1])
    mapping = {"1" : "Genome"}
    annotation = genome.read(open(sys.argv[2]), format="broad", mapping=mapping)

    print ("Chromosome\tStart\tStop\tHeight\tUpstream\tDownstream\tGenic")
    for p in poisson.parse(peaks_handle):
        reg_genes = p.get_regulated_genes(annotation)
        [genic, upstream, downstream] = split_features(reg_genes)

        print("%s\t%d\t%d\t%d\t%s\t%s\t%s" % (p.chromosome, p.chrom_start,
                                              p.chrom_end, p.score(),
                                              ",".join(upstream),
                                              ",".join(downstream),
                                              ",".join(genic)))

