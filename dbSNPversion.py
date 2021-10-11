import pandas as pd
import requests, sys, json, psycopg2, re

# check if RSIDs are valid for build 155 of dbSNP

# Connect to local postgres dbsnp DB.
def connect_dbSNP():
    conn = psycopg2.connect("dbname=dbsnp user=postgres")
    cur = conn.cursor()

    return cur


def dbSNPversion(cur, id):
    query = "SELECT * FROM mergeitems WHERE old_snp_id='" + id + "';"
    cur.execute(query)
    records = cur.fetchall()

    if records:
        return records
    else:
        return -1


def checkDataset():
    # currently using BM, just replace
    file = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/output/biomarkers.tsv', sep='\t', header=0)
    cur = connect_dbSNP()

    for row in file.itertuples():
        result = dbSNPversion(cur, str(row[2]))
        if result != "-1" and result != -1:
            print(str(row[2]) + "= " + str(result))


checkDataset()
