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

    def get_chip_interactions(self, operons = True, min_score = 12.0, 
                             lsr2_threshold = 2.0):
        """Gets the ChIP interactions from the database.  This uses one of
        the two views Anna set up, either chip_interactions_bd_operons if
        the `operons` parameter is true, or chip_interactions_bd otherwise.
        
        Peaks can be filtered by either lsr2 score (default 2.0 or greater)
        and also by height over median (default 12.0)
        """
        # Creating a new array, as we don't need to store everything
        results = []
        
        cursor = self.cursor(MySQLdb.cursors.DictCursor)
        
        if operons:
            db = "chip_interactions_bd_operons"
        else:
            db = "chip_interactions_db"
        
        cursor.execute("SELECT * from %s WHERE `height/median` > %s \
AND (lsr2_ratio = NULL OR lsr2_ratio > %s)" % (db, min_score, lsr2_threshold))
        
        for r in cursor.fetchall():
            results.append((r["from_gene_locus"], r["gene2"]))
        
        return results
        
        