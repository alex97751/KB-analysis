import csv, time, re

import pandas as pd

from getRSID import getRSID, getRSIDfromDB, getRSIDfromDB_using_locusID, connect_dbSNP
from getRefSeq import getRefSeqId
from parseHGVS import setHGVS
from variant import Variant


def main():
    start = time.time()
    # INSERT FUNCTION
    end = time.time()

    print("Duration: " + str(end - start) + "s ")


def progressBar(current, total):
    barLength = 100
    percent = current / total * 100
    arrow = '-' * int(percent / 100 * barLength - 1) + '>'
    spaces = ' ' * (barLength - len(arrow))

    print(' Progress: [%s%s] %.2f %%' % (arrow, spaces, percent), end='\r')


# samples
def parse_cosmic_DB():
    chunksize = 10 ** 3
    file_input = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/KB/cosmic/CosmicMutantExport.tsv', sep='\t', header=0, chunksize=chunksize)

    file_output = open("../output/cosmic_rsid_DB.tsv", "w")
    file_output.write("# \t Found \t Assembly \t HGVS \t Validity \t Error type \t Mutation \t ENS \t RefSeq \t RSID \n")

    count = 0
    found = 0

    dbSNP = connect_dbSNP()  # connect to local dbSNP postgres DB
    for chunk in file_input:
        for row in chunk.itertuples():
            count += 1

            variant = Variant()
            variant.hgvs = str(row[39])
            variant.assembly = 38

            if True:
            # if row[28] == "y":
                setHGVS(variant)
                getRefSeqId(variant)
                for r in variant.refseq:
                    variant.rsid = getRSIDfromDB(dbSNP, r, variant.posedit)
                    if(variant.rsid != -1):
                        found += 1
                        break

            output = (str(count) + "\t" + str(found) + "\t " + str(variant.assembly) + "\t" + str(variant.hgvs) + "\t" + str(variant.valid) + "\t" + str(variant.hgvs_error) + "\t" + str(variant.posedit) + "\t" + str(variant.ac) + "\t" + str(variant.refseq) + "\t" + str(variant.rsid) + "\n")
            file_output.write(output)

            # print(output)
            progressBar(count, chunksize)
        # break  # DEBUG stop after first chunk


# FINISHED
def parse_oncokb_DB(): # uses ENTREZ ID not refseq or ENS
    file_input = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/KB/oncokb/allAnnotatedVariants.txt', sep='\t', header=0, error_bad_lines=False)

    file_output = open("../output/oncokb_rsid_DB.tsv", "w")
    file_output.write("# \t Found \t Mutation \t Entrez ID \t RSID \n")

    count = 0
    found = 0

    dbSNP = connect_dbSNP()  # connect to local dbSNP postgres DB
    for row in file_input.itertuples():
        count += 1

        variant = Variant()
        variant.entrezid = int(row[3])
        variant.posedit = str(row[6])  # always protein change
        variant.proteinOneLetter()

        variant.rsid = getRSIDfromDB_using_locusID(dbSNP, str(variant.entrezid), variant.posedit)

        if(variant.rsid != -1):
            found += 1

        output = str(count) + "\t" + str(found) + "\t" + str(variant.posedit) + "\t" + str(variant.entrezid) + "\t" + str(variant.rsid) + "\n"

        file_output.write(output)

        # print(output)
        progressBar(count, len(file_input.index))


# FINISHED
def parse_civic_DB():
    file_input = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/KB/civic/nightly-VariantSummaries.tsv', sep='\t', header=0, error_bad_lines=False)

    file_output = open("../output/civic_rsid_DB.tsv", "w")
    file_output.write("# \t Found \t Assembly \t HGVS \t Validity \t Error type \t Mutation \t Transcript IDs \t RSID \n")

    count = 0
    found = 0

    dbSNP = connect_dbSNP()  # connect to local dbSNP postgres DB
    for row in file_input.itertuples():
        count += 1

        variant = Variant()
        hgvs_list = str(row[21]).split(",")

        assembly = re.findall(r'\d+', str(row[15]))
        if len(assembly) > 0:
            variant.assembly = assembly[0]

        for v in hgvs_list:
            variant.hgvs = v
            setHGVS(variant)

            if(True):
                getRefSeqId(variant)
                for r in variant.refseq:
                    variant.rsid = getRSIDfromDB(dbSNP, r, variant.posedit)
                    if(variant.rsid != -1):
                        break

            if(variant.rsid != -1):
                found += 1
                break

        output = str(count) + "\t" + str(found) + "\t " + str(variant.assembly) + "\t" + str(variant.hgvs) + "\t" + str(variant.valid) + "\t" + str(variant.hgvs_error) + "\t" + str(variant.posedit) + "\t" + str(variant.refseq) + "\t" + str(variant.rsid) + "\n"
        file_output.write(output)

        # print(output)
        # print("%0.2f" % (count/len(file_input.index) * 100) + " ", end='', flush=True)
        progressBar(count, len(file_input.index))


# ANALYSE HGVS IN CIVIC
def parse_civic_HGVS():
    file_input = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/KB/civic/nightly-VariantSummaries.tsv', sep='\t', header=0, error_bad_lines=False)

    file_output = open("../output/civic_HGVS.tsv", "w")
    file_output.write("Total \t Variant \t # \t Found Valid \t HGVS \t Found RefSeqs \t RefSeq \t Validity \t Error type \n")

    variants = 0
    totalHGVSstrings = 0
    validHGVS = 0
    totalRefSeqs = 0
    count = 0
    refSeq = False

    for row in file_input.itertuples():
        variants += 1

        variant = Variant()
        hgvs_list = str(row[21]).split(",")

        listcounter = 0
        for v in hgvs_list:
            count += 1
            listcounter += 1
            totalHGVSstrings += 1
            variant.hgvs = v
            setHGVS(variant)
            if(variant.valid):
                validHGVS += 1
            if(str(variant.hgvs).startswith("N")):
                totalRefSeqs += 1
                refSeq = True
            else:
                refSeq = False


            output = str(count) + "\t" + str(variants) + "\t" + str(listcounter) + "\t" + str(validHGVS) + "\t " + str(variant.hgvs) + "\t" + str(totalRefSeqs) + "\t" + str(refSeq) + "\t" + str(variant.valid) + "\t" + str(variant.hgvs_error) + "\n"
            file_output.write(output)

        # print(output)
        # print("%0.2f" % (count/len(file_input.index) * 100) + " ", end='', flush=True)
        progressBar(variants, len(file_input.index))

    print("Variants: " + str(variants))
    print("total HGVS Strings: " + str(totalHGVSstrings))
    print("valid RefSeq Strings: " + str(totalRefSeqs))
    print("valid HGVS Strings: " + str(validHGVS))


# FINISHED
def parse_biomarkers_DB():
    file_input = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/KB/biomarkers/cgi_biomarkers_latest/cgi_biomarkers_per_variant.tsv', sep='\t', header=0)

    file_output = open("../output/biomarkers_rsid_DB.tsv", "w")
    file_output.write("# \t Found \t Assembly \t HGVS \t Validity \t Error type \t Mutation \t ENS \t RefSeq \t RSID \n")

    count = 0
    found = 0

    dbSNP = connect_dbSNP()  # connect to local dbSNP postgres DB
    for row in file_input.itertuples():
        count += 1

        variant = Variant()
        variant.hgvs = str(row[27]) + ":" + str(row[21])
        variant.assembly = 38

        setHGVS(variant)

        if(True):
            getRefSeqId(variant)
            for r in variant.refseq:
                variant.rsid = getRSIDfromDB(dbSNP, r, variant.posedit)
                if(variant.rsid != -1):
                    found += 1
                    break

        output = str(count) + "\t" + str(found) + "\t " + str(variant.assembly) + "\t" + str(variant.hgvs) + "\t" + str(variant.valid) + "\t" + str(variant.hgvs_error) + "\t" + str(variant.posedit) + "\t" + str(variant.ac) + "\t" + str(variant.refseq) + "\t" + str(variant.rsid) + "\n"

        file_output.write(output)

        #print("%0.2f" % (count/len(file_input.index) * 100) + " ", end='', flush=True)
        progressBar(count, len(file_input.index))


# FINISHED
def parse_docm_DB():
    file_input = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/KB/docm/variants.tsv', sep='\t', header=0)

    file_output = open("../output/docm_rsid_DB.tsv", "w")
    file_output.write("# \t Found \t Assembly \t HGVS \t Validity \t Error type \t Mutation \t ENS \t RefSeq \t RSID \n")

    count = 0
    found = 0

    dbSNP = connect_dbSNP()  # connect to local dbSNP postgres DB
    for row in file_input.itertuples():
        count += 1

        variant = Variant()
        variant.hgvs = str(row[1])

        assembly = re.findall(r'\d+', str(row[7]))  # assembly used (always 38 here!!)
        if len(assembly) > 0:
            variant.assembly = assembly[0]

        setHGVS(variant)

        if(True):
            getRefSeqId(variant)
            for r in variant.refseq:
                variant.rsid = getRSIDfromDB(dbSNP, r, variant.posedit)
                if(variant.rsid != -1):
                    found += 1
                    break

        output = str(count) + "\t" + str(found) + "\t " + str(variant.assembly) + "\t" + str(variant.hgvs) + "\t" + str(variant.valid) + "\t" + str(variant.hgvs_error) + "\t" + str(variant.posedit) + "\t" + str(variant.ac) + "\t" + str(variant.refseq) + "\t" + str(variant.rsid) + "\n"

        file_output.write(output)

        # print("%0.2f" % (count/len(file_input.index) * 100) + " ", end='', flush=True)
        progressBar(count, len(file_input.index))

    file_output.close()


# sample
def parse_dbNSFP4():
    file_input = open("../variant_dbs_small/dbNSFP4.1a_variant.chr1.tsv")
    read_tsv = csv.reader(file_input, delimiter="\t")
    for row in read_tsv:
        variant = Variant()
        variant.ac = row[14]
        variant.posedit = row[22]
        # variant.ac = getRefSeqId(variant.ac)
        # variant.rsid = getRSIDfromDB(str(variant.ac), str(variant.posedit))
        variant.rsid = getRSID(str(variant.ac) + ":" + str(variant.posedit))
        print(str(variant.ac) + ":" + str(variant.posedit) + " - " + str(variant.rsid))

# sample
def parse_clinvar_DB():
    file_input = open("../variant_dbs_small/clinvar_variant_summary.txt")
    read_tsv = csv.reader(file_input, delimiter="\t")
    count = 0
    foundHgvs = 0
    found = 0

    for row in read_tsv:
        count += 1
        variant = Variant()
        variant.hgvs = row[2]
        r = setHGVS(variant.hgvs)

        if(r != -1):
            variant.rsid = getRSIDfromDB(r[0], r[1])

            foundHgvs += 1
            if(variant.rsid != "-1"):
                found += 1
            print(str(count) + "( " + str(foundHgvs) + "/" + str(found) + "): " + str(variant.hgvs) + " - " + str(variant.rsid))

    file_input.close()


# FINISHED
def parse_docm_API():
    file_input = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/KB/docm/variants.tsv', sep='\t', header=0)
    count = 0
    found = 0

    file_output = open("../output/docm_rsid_API.tsv", "w")
    file_output.write("# \t Found \t Assembly \t HGVS \t RSID \n")

    for row in file_input.itertuples():
        variant = Variant()
        variant.hgvs = str(row[1])
        variant.assembly = str(re.findall(r'\d+', str(row[7]))[0])
        variant.rsid = getRSID(str(row[1]), variant.assembly)

        count += 1
        if(variant.rsid != "-1"): found += 1

        output = str(count) + "\t" + str(found) + "\t " + str(variant.assembly) + "\t" + str(variant.hgvs) + "\t" + str(variant.rsid) + "\n"

        # print(output)
        file_output.write(output)

        progressBar(count, len(file_input.index))

    file_output.close()
    file_input.close()


# FINISHED
def parse_civic_API():
    file_input = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/KB/civic/nightly-VariantSummaries.tsv', sep='\t', header=0, error_bad_lines=False)

    file_output = open("../output/civic_rsid_API.tsv", "w")
    file_output.write("# \t Found \t Assembly \t HGVS \t RSID \n")

    count = 0
    found = 0

    for row in file_input.itertuples():
        count += 1

        variant = Variant()
        hgvs_list = str(row[21]).split(",")

        assembly = re.findall(r'\d+', str(row[15]))
        if len(assembly) > 0:
            variant.assembly = assembly[0]

        for v in hgvs_list:
            variant.hgvs = v
            variant.rsid = getRSID(str(variant.hgvs), variant.assembly)

            if(variant.rsid != "-1"):
                found += 1
                break

        output = str(count) + "\t" + str(found) + "\t " + str(variant.assembly) + "\t" + str(variant.hgvs) + "\t" + str(variant.rsid) + "\n"
        file_output.write(output)

        # print(output)
        # print("%0.2f" % (count/len(file_input.index) * 100) + " ", end='', flush=True)
        progressBar(count, len(file_input.index))

#sample
def parse_cosmic_API():
    chunksize = 10 ** 3
    file_input = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/KB/cosmic/CosmicMutantExport.tsv', sep='\t', header=0, chunksize=chunksize)

    count = 0
    found = 0

    file_output = open("../output/cosmic_rsid_API.tsv", "w")
    file_output.write("# \t Found \t HGVS \t RSID \n")

    for chunk in file_input:
        for row in chunk.itertuples():
            variant = Variant()
            variant.hgvs = str(row[39])
            #print(variant.hgvs)
            variant.rsid = str(getRSID(variant.hgvs, "38"))
            # use hgvsp if hgvsc is not found
            if(variant.rsid == "-1"):
                variant.hgvs = str(row[38])
                #print("NEW: " + variant.hgvs)
                variant.rsid = str(getRSID(variant.hgvs, "38"))

            count += 1
            if(variant.rsid != "-1"): found += 1

            output = (str(count) + "\t" + str(found) + "\t" + str(variant.hgvs) + "\t" + str(variant.rsid) + "\n")
            file_output.write(output)

            progressBar(count, chunksize)
        break

    file_output.close()
    file_input.close()


#sample
def parse_clinvar_API():
    file_input = open("../variant_dbs_small/clinvar_variant_summary.txt")
    read_tsv = csv.reader(file_input, delimiter="\t")
    count = 0
    found = 0
    different_rsids = 0
    variants = list()

    file_output = open("../output/clinvar_rsid_API.txt", "w")

    for row in read_tsv:
        #if count > 50: break
        variant = Variant()
        variant.hgvs = row[2].split(" ")  # TODO neccesary?
        variant.rsid = str(row[9])

        # rsid is not provided for this variant
        if(variant.rsid == "-1"):
            for i in range(0, len(variant.hgvs)):
                variant.rsid = str(getRSID(variant.hgvs[i]))
                if(variant.rsid != "-1"): break;
        else: # compare with RSID in file
            for i in range(0, len(variant.hgvs)):
                compare = str(getRSID(variant.hgvs[i]))
                if(compare != "rs"+variant.rsid): different_rsids += 1
                if(variant.rsid != "-1"): break;

        variants.append(variant)

        count += 1
        if(variant.rsid != "-1"): found += 1
        print(str(count) + "( " + str(found) + "): " + variant.hgvs[0] + " - " + variant.rsid)
        file_output.write(str(count) + "( " + str(found) + "): " + variant.hgvs[0] + " - " + variant.rsid + "\n")

    file_output.close()
    file_input.close()
    print("different RSIDs: " + str(different_rsids) + " / " + str(count))


# FINISHED
def parse_biomarkers_API():
    file_input = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/KB/biomarkers/cgi_biomarkers_latest/cgi_biomarkers_per_variant.tsv', sep='\t', header=0)

    file_output = open("../output/biomarkers_rsid_API.tsv", "w")
    file_output.write("# \t Found \t Assembly \t HGVS \t RSID \n")

    count = 0
    found = 0

    for row in file_input.itertuples():
        count += 1

        variant = Variant()
        variant.hgvs = str(row[27]) + ":" + str(row[21])
        variant.assembly = 38
        variant.rsid = getRSID(str(variant.hgvs), variant.assembly)

        if(variant.rsid != "-1"): found += 1

        output = str(count) + "\t" + str(found) + "\t " + str(variant.assembly) + "\t" + str(variant.hgvs) + "\t" + str(variant.rsid) + "\n"
        file_output.write(output)

        progressBar(count, len(file_input.index))


# FINISHED
def parse_oncokb_API(): # uses ENTREZ ID not refseq or ENS
    file_input = pd.read_csv('/Users/alexharrisson/Documents/Uni/BA/KB/oncokb/allAnnotatedVariants.txt', sep='\t', header=0, error_bad_lines=False)

    file_output = open("../output/oncokb_rsid_API.tsv", "w")
    file_output.write("# \t Found \t Mutation \t HGVS ID \t RSID \n")

    count = 0
    found = 0

    for row in file_input.itertuples():
        count += 1

        variant = Variant()
        variant.ac = str(row[1])
        variant.posedit = str(row[6])  # always protein change
        variant.proteinOneLetter()
        variant.hgvs = str(variant.ac) + ":" + str(variant.posedit)

        variant.rsid = getRSID(str(variant.hgvs), variant.assembly)
        if(variant.rsid != "-1"):
            found += 1

        output = str(count) + "\t" + str(found) + "\t" + str(variant.posedit) + "\t" + str(variant.hgvs) + "\t" + str(variant.rsid) + "\n"

        file_output.write(output)

        # print(output)
        progressBar(count, len(file_input.index))


if __name__ == "__main__":
    main()
