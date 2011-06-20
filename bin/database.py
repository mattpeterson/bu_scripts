#!/usr/bin/env python

import MySQLdb
from collections import defaultdict

class DatabaseConnection(object):
    """
    A connection to the TB network database
    """
    def __init__(self, host = "tuberculosis.bu.edu", db  = "tbnetwork_v3",
                 user = "petersmw", password = "MP_bnl8182"):
        self._connection = MySQLdb.connect(host = host, user = user,
                                           passwd = password, db = db)
    
    def cursor(self):
        """
        Returns the cursor for writing custom queries.
        """
        return self._connection.cursor()
    
    def get_orthogroups(self, org_code = "h37r"):
        """
        Gets the orthogroups for an organism
        
        :param org_code: The organism code to return orthogroups
        :type org_code: string
        """
        orthogroups = defaultdict(list)
        
        cursor = self.cursor()
        
        results = cursor.execute("SELECT g.locus, g.orthogroups FROM genes \
as g WHERE g.org_code LIKE %s", (org_code,))
        
        for result in results:
            orthogroups[results[0]].append(results[1])
        
        return orthogroups
    