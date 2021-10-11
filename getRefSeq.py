import mysql.connector

# connect to homo sapiens DB 37 or 38 to get a valid RefSeq ID

# Connect to local mysql DB for ENS->RefSeq mapping.
def connect(variant):
    if variant.assembly == 37:
        variant.mydb = mysql.connector.connect(host="localhost", user="root", password="", database="homo_sapiens_core_75_37")
        variant.mycursor = variant.mydb.cursor()
    else:  # use assembly 38
        variant.mydb = mysql.connector.connect(host="localhost", user="root", password="", database="homo_sapiens_core_104_38")
        variant.mycursor = variant.mydb.cursor()

    return


# Query mysql DB to map ENS->RefSeq.
def getRefSeqId(variant):
    if(variant.ac is not None):
        connect(variant)
        ensembl_id = variant.ac.split(".")[0]
        query = "SELECT xref.display_label FROM transcript, object_xref, xref,external_db WHERE transcript.transcript_id = object_xref.ensembl_id AND object_xref.ensembl_object_type = 'Transcript' AND object_xref.xref_id = xref.xref_id AND xref.external_db_id = external_db.external_db_id AND external_db.db_name = 'RefSeq_mRNA' AND transcript.stable_id ='" + ensembl_id + "';"
        variant.mycursor.execute(query)
        myresult = variant.mycursor.fetchall()

        for refseq in myresult:
            variant.refseq.append(refseq[0])

    return

# DEBUG
# ensembl_id = "ENST00000497140"
# getRefSeqId(ensembl_id)
