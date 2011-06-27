#!/usr/bin/env python

import MySQLdb
import MySQLdb.cursors
from collections import defaultdict

class DatabaseConnection(object):
    """
    A connection to the TB network database
    """
    def __init__(self, host = "tuberculosis.bu.edu", db  = "tbnetwork_v3",
                 user = "petersmw", password = "MP_bnl8182"):
        self._connection = MySQLdb.connect(host = host, user = user,
                                           passwd = password, db = db)
    
    def cursor(self, cursorclass = MySQLdb.cursors.Cursor):
        """
        Returns the cursor for writing custom queries.
        """
        return self._connection.cursor(cursorclass = cursorclass)
    
    def get_orthogroups(self, org_code = "h37r", key = "gene"):
        """
        Gets the orthogroups for an organism
        
        :param org_code: The organism code to return orthogroups
        :type org_code: string
        """
        orthogroups = defaultdict(list)
        
        cursor = self.cursor()
        
        results = cursor.execute("SELECT g.locus, g.orthogroups FROM genes \
as g WHERE g.org_code LIKE %s", (org_code,))
        
        if key == "gene":
            for result in results:
                orthogroups[results[0]].append(results[1])
        elif key == "orthogroup":
            for result in results:
                orthogroups[results[1]].append(results[0])
        else:
            raise ValueError("Key must be 'gene' or 'orthogroup'.")
        
        return orthogroups

    def get_chip_interactions(self, operons = True, min_score = 4.0, 
                             lsr2_threshold = 2.0):
        """Gets the ChIP interactions from the database.  This uses one of
        the two views Anna set up, either chip_interactions_bd_operons if
        the `operons` parameter is true, or chip_interactions_bd otherwise.
        
        Peaks can be filtered by either lsr2 score (default 2.0 or greater)
        and also by height over median (default 12.0)
        """
        # Creating a new array, as we don't need to store everything
        results = []
        
        cursor = self.cursor()
        
        # Using the chip_interactions_simple view.
        #cursor.execute("SELECT from_gene_locus, to_gene_locus from \
#chip_interactions_simple")
        cursor.execute("SELECT from_gene_locus, to_gene_locus from \
chip_interactions_bd WHERE `height/median` > 5 AND (lsr2_ratio IS \
NULL or lsr2_ratio > 2)")
        for r in cursor.fetchall():
            results.append((r[0], r[1]))
        
        return results
    
    
    def get_clr_scores(self, gene = None, ignore_zeros = False):
        """Gets the CLR scores for all pairs of genes.  Returns a dictionary,
        keyed on tuples.  To save space, don't return anything w/ 0.0
        """
        results = {}
        
        cursor = self.cursor()
        
        if ignore_zeros:
            cursor.execute("SELECT * from oe_clr WHERE z_score > 0 \
AND from_gene_locus = %s", (gene,))
        else:
            cursor.execute("SELECT * from oe_clr WHERE from_gene_locus = %s",
                           (gene,))
        
        for r in cursor.fetchall():
            results[(r[0], r[1])] = r[2]
            
        return results
    
    def get_lit_interactions(self, gene):
        """Retrieve literature curated interactions from the database"""
        interactions = set()
        
        cursor = self.cursor()
        
        cursor.execute("SELECT from_gene_locus, to_gene_locus FROM \
interactions WHERE type = \"Curation\" AND from_gene_locus = %s", (gene,))
        
        # Store everything as a tuple.  Easier that way.
        for r in cursor.fetchall():
            interactions.add( (r[0], r[1]) )
        
        return interactions