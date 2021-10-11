import pandas as pd

# Compare Results from the API and the manual search approach

# 1 db
def compareBM():
    file_DB = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/output/biomarkers_rsid_DB.tsv', sep='\t', header=0)
    file_API = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/output/biomarkers_rsid_API.tsv', sep='\t', header=0)

    list_db = list()
    list_api = list()
    for row in file_DB.itertuples():
        list_db.append(row[10])

    for row in file_API.itertuples():
        list_api.append(row[5])

    for i in range(len(list_db)):
        if(list_db[i] != list_api[i]):
            print(str(list_db[i]) + " - " + str(list_api[i]))
        i += 1


# 3 db
# 1 missmatch
def compareCivic():
    file_DB = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/output/civic_rsid_DB.tsv', sep='\t', header=0)
    file_API = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/output/civic_rsid_API.tsv', sep='\t', header=0)

    list_db = list()
    list_api = list()
    for row in file_DB.itertuples():
        list_db.append((row[9], row[4]))

    for row in file_API.itertuples():
        list_api.append((row[5], row[4]))

    for i in range(len(list_db)):
        if(list_db[i][0] != list_api[i][0]):
            if(list_db[i][0] > 0 and list_api[i][0] > 0):
                print(str(list_db[i]) + " - " + str(list_api[i]))
        i += 1


# 1 missmatch
def compareDocm():
    file_DB = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/output/docm_rsid_DB.tsv', sep='\t', header=0)
    file_API = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/output/docm_rsid_API.tsv', sep='\t', header=0)

    list_db = list()
    list_api = list()
    for row in file_DB.itertuples():
        list_db.append((row[10], row[4]))

    for row in file_API.itertuples():
        list_api.append((row[10], row[4]))

    for i in range(len(list_db)):
        if(list_db[i][0] != list_api[i][0]):
            if(list_db[i][0] > 0 and list_api[i][0] > 0):
                print(str(list_db[i]) + " - " + str(list_api[i]))
        i += 1


# 357 only in db
# 26 missmatch
def compareOncoKB():
    file_DB = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/output/oncokb_rsid_DB.tsv', sep='\t', header=0)
    file_API = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/output/oncokb_rsid_API.tsv', sep='\t', header=0)

    list_db = list()
    list_api = list()
    for row in file_DB.itertuples():
        list_db.append([row[5], row[4]])

    for row in file_API.itertuples():
        list_api.append([row[5], row[4]])

    count = 0
    for i in range(len(list_db)):
        if(list_db[i] != list_api[i]):  # not similar entries
            if(list_db[i][0] > 0 and list_api[i][0] > 0):  # find mismatches
                count += 1
                print(str(list_db[i][0]) + " - " + str(list_api[i][0]))
        i += 1
    print(count)


# compare API and DB results for each KB
compareCivic()
