#!/usr/bin/env/python
"""
A script to read in the results of a BCG ChIP experiment, and calculate 
enrichment using the MTB COG categories
"""
import sys
from database import DatabaseConnection

def main():
    # First, let's get BCG orthogroups (gene keyed)
    db_connection = DatabaseConnection()
    bcg_to_orthogroup = db_connection.get_orthogroups("bbcg")
    # Now, the TB orthogroups
    orthogroup_to_mtb = db_connection.get_orthogroups("h37r", "orthogroups")
    
    results_file = sys.argv[1]    

if __name__ == "__main__":
    main()

