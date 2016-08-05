import sqlite3
import ete2
# sys.argv[1] = $the ncbi.db database
# the dabase is required to parse the taxonomy from a gi_id



for line in open(sys.argv[2], 'r'):
    data = line.split()
    
