#!/usr/bin/env python
"""
Pulls the curated results from the database and compares them to the 
ChIP-seq interactions
"""
from database import DatabaseConnection
import sys

def get_matching_pct(curated_interactions, chip_interactions):
    total_curated = len(curated_interactions)
    
    hits = 0
    for i in curated_interactions:
        if i in chip_interactions:
            hits += 1
    
    return float(hits) / total_curated

def read_balazi_network(file_handle):
    interactions = set()
    
    for line in file_handle:
        tokens = line.rstrip("\r\n").split("\t")
        interactions.add( (tokens[2], tokens[5]))
    
    return interactions

def main():
    connection = DatabaseConnection()
    
    balazi_interactions = read_balazi_network(open(sys.argv[1]))
    chip_interactions = connection.get_chip_interactions()
    from_genes = set([ i[0] for i in chip_interactions ])
    
    for g in from_genes:
        curated_interactions = filter(lambda x: x[0] == g, balazi_interactions)
        
        if len(curated_interactions) == 0:
            print("%s: None Curated" % g)
            continue
        
        print("%s: %f" % (g, get_matching_pct(curated_interactions, 
                                              chip_interactions)))
        

if __name__ == "__main__":
    main()
