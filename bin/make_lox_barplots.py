#!/usr/bin/env python
"""
Makes Barplots of ChIP-seq and RNA-seq values.  Used as a one-off script for
Nicole, but could be put into a more functional, general script.

Matthew Peterson
June 20, 2011
"""
import sys
import transnet.transcriptomics.lox as lox
from matplotlib import pyplot as plt
import numpy as np

def read_chip(chip_handle):
    
    values = [{} for i in range(3)]
    errors = [{} for i in range(3)]
    
    for line in chip_handle:
        tokens = line.rstrip("\r\n").split("\t")
        print tokens
        values[0][tokens[0]] = float(tokens[2])
        errors[0][tokens[0]] = (float(tokens[3]), float(tokens[4]))
        
        values[1][tokens[0]] = float(tokens[5])
        errors[1][tokens[0]] = (float(tokens[6]), float(tokens[7]))
        
        values[2][tokens[0]] = float(tokens[8])
        errors[2][tokens[0]] = (float(tokens[9]), float(tokens[10]))
    
    return (values, errors)

def make_barplot(axis, experiments, locus, plot_groups, names):
    """
    Produces a bar plot for the set of experiments, grouped by the plot_groups
    option.
    """
    # Steal these from the Matplotlib example
    # (http://matplotlib.sourceforge.net/examples/api/barchart_demo.html)
    width = 0.35
    
    rects = []
    colors = ["b", "r"]
    
    for i, group in enumerate(plot_groups):
        values = []
        up_errors = []
        down_errors = []
        
        indices = np.arange(len(group)) + i * width
        for v in group:
            values.append(experiments[v][locus].expression_level)
            up_errors.append(experiments[v][locus].upper_confidence)
            down_errors.append(experiments[v][locus].lower_confidence)
        
        errors = [up_errors, down_errors]
    
        print values
        print indices
        
        rects.append(plt.bar(indices, np.array(values), width, yerr=errors, 
                             color=colors[i]))
        
    axis.set_xticks(np.arange(len(plot_groups[1])) + width)
    axis.set_xticklabels(names)
    axis.set_ylabel("Normalized Expression")
    
    axis.legend([rect[0] for rect in rects], ["WT", "dCSP"])
    

def make_chip_barplot(axis, values, errors, locus):
    width = 0.7
    indices = [0,1,2]
    labels = ["Dark", "15m", "60m"]
    
    plot_values = []
    up_errors = []
    down_errors = []
    
    for i in range(3):
        plot_values.append(values[i][locus])
        up_errors.append(errors[i][locus][0])
        down_errors.append(errors[i][locus][1])
    
    plot_errors = [up_errors, down_errors]
    
    plt.bar(indices, plot_values, width, yerr=plot_errors)
    
    axis.set_xticks(np.array(indices) + 0.35)
    axis.set_xticklabels(labels)
    axis.set_ylabel('ChIP tags (RPKM)')
    
def main():

    # The groups in which everything should be plotted
    plot_groups = [ ["WT-dark", "WT-15m", "WT-60m"],
                  ["dCSP-dark", "dCSP-15m", "dCSP-60m"] ]
    names = ["Dark", "15m", "60m"]

    # Read in the LOX experiment
    rna_experiments = lox.read(open(sys.argv[1]))
    
    # Tons 'o kludge to get the ChIP-seq stuff plotted w/ scale
    (values, errors) = read_chip(open(sys.argv[2]))

    with open(sys.argv[3]) as gene_list_handle:
        for line in gene_list_handle:
            (locus, description) = line.rstrip("\r\n").split("\t")
            
            fig = plt.figure()
            ax1 = fig.add_subplot(211)
            make_barplot(ax1, rna_experiments, locus, plot_groups, names)
            
            ax2 = fig.add_subplot(212)
            make_chip_barplot(ax2, values, errors, locus)
            
            ax1.set_title("%s: %s" % (locus, description))
            
            fig.savefig("./plots/%s.png" % locus, dpi=150)
            
            plt.close()

if __name__ == "__main__":
    main()
