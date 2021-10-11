import pandas as pd

# Compare RSIDs provided by CIViC with RSIDs the manual search and API found

def main():
    file_input = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/output/civic_rsid_db.tsv', sep='\t', header=0, error_bad_lines=False)
    rsids = comapreCivic()
    count = 0
    same = 0
    missmatch = 0

    api = 5
    db = 9
    index = db

    print("Found with DB: 332")
    print("Found with API: 487")
    print("Provided by CIViC: " + str(len(rsids)))

    for row in file_input.itertuples():
        if str(row[index]) == str(rsids[count]):
            if int(row[index]) > 0 and int(rsids[count]) > 0:
                same += 1

        if str(row[index]) != str(rsids[count]):
            if int(row[index]) > 0 and int(rsids[count]) > 0:
                print(row[index], rsids[count])
                missmatch += 1

        count += 1

    print("Identical: " + str(same))
    print("Missmatch: " + str(missmatch))


def comapreCivic():
    file_input = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/KB/civic/nightly-VariantSummaries.tsv', sep='\t', header=0, error_bad_lines=False)

    count = 0
    found = 0
    rsids = list()

    for row in file_input.itertuples():
        rsids.append("-1")
        alias_list = str(row[26]).split(",")

        for v in alias_list:
            if v.startswith("RS"):
                found += 1
                rsids[count] = v.replace("RS", "")

        count += 1

    print("Found RSIDs in CIViC: " + str(found))
    return rsids


main()
