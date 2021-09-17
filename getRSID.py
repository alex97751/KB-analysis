import requests, sys, json, psycopg2, re


# Return RSID of a hgvs variant from ensamble REST api.
def getRSID(hgvs, assembly):
    if(assembly == "37"):
        server = "https://grch37.rest.ensembl.org"
    else:
        server = "https://rest.ensembl.org"

    ext = "/variant_recoder/human/" + hgvs + "?"

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        return "-1"  # input error
    else:
        decoded = r.json()
        try:
            #print(decoded)
            type = str(decoded)[3] # ugly work around
            rsid = decoded[0][type]['id'][0]  # assuming rsid is at first position
            if(rsid[0] == "r" and rsid[1] == "s"):
                # print(rsid)
                return str(re.findall(r'\d+', rsid)[0])
            else:
                return "-1"  # not rsid, probably COSV
        except KeyboardInterrupt:
            sys.exit(0)
        except Exception:
            return "-1"


# Connect to local postgres dbsnp DB.
def connect_dbSNP():
    conn = psycopg2.connect("dbname=dbsnp user=postgres")
    cur = conn.cursor()

    return cur


# Return RSID of a hgvs variant from local dbsnp DB.
def getRSIDfromDB(cur, ac, posedit):
    ac_without_version = ac.split(".")[0]  # TODO ??
    query = "SELECT snp_id FROM hgvs WHERE hgvs='" + posedit + "' AND refseq LIKE'" + ac_without_version + "%';"
    cur.execute(query)
    records = cur.fetchall()

    if records:
        return records[0][0]  # return first found rsid
    else:
        return -1


# Return RSID of a hgvs variant from local dbsnp DB using lous id not transcript.
def getRSIDfromDB_using_locusID(cur, entrezid, posedit):
    query = "SELECT snp_id FROM hgvs WHERE hgvs='" + posedit + "' AND locus_id='" + entrezid + "';"
    cur.execute(query)
    records = cur.fetchall()

    if records:
        return records[0][0]  # return first found rsid
    else:
        return -1

# DEBUG
# getRSID("ENST00000490903.5c.1402+1908G>A")
# getRSIDfromDB("NM_001386501.1", "c.6643T>A")
