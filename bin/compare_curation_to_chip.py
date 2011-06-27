#!/usr/bin/env python
"""
Pulls the curated results from the database and compares them to the 
ChIP-seq interactions
"""
from database import DatabaseConnection
import sys

def get_matching_pct(curated_interactions, chip_interactions):
    """
    Returns percentage of matching
    """
    total_curated = len(curated_interactions)
    
    hits = 0
    for i in curated_interactions:
        if i in chip_interactions:
            hits += 1
    
    return float(hits)

def read_balazi_network(file_handle):
    interactions = set()
    
    for line in file_handle:
        tokens = line.rstrip("\r\n").split("\t")
        if "1" in tokens[4] or "2" in tokens[4]:
            interactions.add( (tokens[1], tokens[3]))
    
    return interactions

def main():
    """
    Main logic
    """
    connection = DatabaseConnection()
    
    balazi_interactions = read_balazi_network(open(sys.argv[1]))

    chip_interactions = connection.get_chip_interactions()
    from_genes = set([ i[0] for i in chip_interactions ])
    
    for g in from_genes:
        curated_interactions = filter(lambda x: x[0] == g, balazi_interactions)
        chip_interactions_g = filter(lambda x: x[0] == g, chip_interactions)
        
        if len(curated_interactions) == 0:
            print("%s: None Curated" % g)
            continue
            
        matches = get_matching_pct(curated_interactions, chip_interactions)
        
        print("%s: %d/%d (%f)" % (g, matches, len(curated_interactions),
                                              float(matches)/len(curated_interactions)))
        

if __name__ == "__main__":
    main()
