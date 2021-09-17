import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles


def main():
    # variants in each kb
    biomarkers = 1442
    civic = 3318
    oncokb = 4915
    docm = 1365
    clinvar = 1867779

    path = 'Users/aharrisson/Documents/uni/BA/output/copy'

    # DOCM
    # file = pd.read_csv('/Users/alexharrisson/Documents/uni/BA/output/docm_rsid_DB.tsv', dtype='unicode', sep='\t', header=0)
    # hashMapDocm_DB = getDict(file, 10)
    # file = pd.read_csv('/Users/alexharrisson/Documents/uni/BA/output/docm_rsid_API.tsv', dtype='unicode', sep='\t', header=0)
    # hashMapDocm_API = getDict(file, 10)

    # Biomarkers
    # file = pd.read_csv('/Users/alexharrisson/Documents/uni/BA/output/biomarkers_rsid_DB.tsv', dtype='unicode', sep='\t', header=0)
    # hashMapBM_DB = getDict(file, 10)
    # file = pd.read_csv('/Users/alexharrisson/Documents/uni/BA/output/biomarkers_rsid_API.tsv', dtype='unicode', sep='\t', header=0)
    # hashMapBM_API = getDict(file, 5)

    # ClinVar
    # file = pd.read_csv('/Users/alexharrisson/Documents/uni/BA/output/clinvar_rsid.tsv', dtype='unicode', sep='\t', header=0)
    # hashMapClinVar = getDict(file, 3)

    # CiVic
    # file = pd.read_csv('/Users/alexharrisson/Documents/uni/BA/output/civic_rsid_DB.tsv', dtype='unicode', sep='\t', header=0)
    # hashMapCivic_DB = getDict(file, 9)
    # file = pd.read_csv('/Users/alexharrisson/Documents/uni/BA/output/civic_rsid_API.tsv', dtype='unicode', sep='\t', header=0)
    # hashMapCivic_API = getDict(file, 5)

    # # OncoKB
    # file = pd.read_csv('/Users/alexharrisson/Documents/uni/BA/output/oncokb_rsid_DB.tsv', dtype='unicode', sep='\t', header=0)
    # hashMapOncoKb_DB = getDict(file, 5)
    # file = pd.read_csv('/Users/alexharrisson/Documents/uni/BA/output/oncokb_rsid_API.tsv', dtype='unicode', sep='\t', header=0)
    # HashMapOncoKb_API = getDict(file, 5)

    # Cosmic
    # file = pd.read_csv('/Users/alexharrisson/Desktop/cosmic_rsid_DB.tsv', dtype='unicode', sep='\t', header=0)
    # hashMapCosmic_DB = getDict(file, 8)
    # file = pd.read_csv('/Users/alexharrisson/Desktop/cosmic_rsid_API.tsv', dtype='unicode', sep='\t', header=0)
    # hashMapCosmic_API = getDict(file, 4)

    # dicts = [hashMapClinVar] # combine these into one file

    # checkSimilarity(dicts)
    # combineBoth(dicts, "clinvar")
    # hashMapOverlap = getOverlapDict(dicts)
    #overlap = getOverlapInt(hashMapOverlap, len(dicts))

    listOfKBs = ["docm", "oncokb", "biomarkers", "civic"]
    lenOfKBs = [docm, oncokb, biomarkers, civic]
    overlap = getOverlapDict(listOfKBs)
    overlapInt = getOverlapInt(overlap, len(listOfKBs))

    print(str(overlapInt) + "/" + str(len(overlap)) + " => " + str(overlapInt / (len(overlap)) * 100) + "%")

    # printVenn2(listOfKBs, lenOfKBs, overlapInt)
    # printVenn2(listOfKBs, lenOfKBs, overlapInt)


def printVenn2(listOfKBs, lenOfKBs, overlapInt):
    plt.figure(figsize=(5, 5))
    venn2(subsets=(lenOfKBs[0], lenOfKBs[1], overlapInt), set_labels=(listOfKBs[0], listOfKBs[1]))
    venn2_circles(subsets=(lenOfKBs[0], lenOfKBs[1], overlapInt), linewidth=1, color='k')

    #plt.show()
    plt.savefig('../thesis/graphs/overlap/' + listOfKBs[0] + "_" + listOfKBs[1]+'.png')


def printVenn3(listOfKBs, lenOfKBs, overlapInt):
    plt.figure(figsize=(5, 5))
    venn3(subsets=(lenOfKBs[0], lenOfKBs[1],lenOfKBs[2],1,2,3, overlapInt), set_labels=(listOfKBs[0], listOfKBs[1], listOfKBs[2]))
    venn3_circles(subsets=(lenOfKBs[0], lenOfKBs[1],lenOfKBs[2], 1, 2, 3, overlapInt), linewidth=1, color='k')

    plt.show()
    # plt.savefig('../thesis/graphs/overlap/' + listOfKBs[0] + "_" + listOfKBs[1]+'.png')


# count overlap between KBs
def getOverlapInt(hashMapOverlap, number):
    overlap = 0
    for i in hashMapOverlap.items():
        if(i[1] == number):
            print(i[0])
            overlap += 1

    return overlap


# combine all KB dict into one dict (hashMap)
def getOverlapDict(dicts):
    hashMapOverlap = {}
    for i in dicts:
        file = pd.read_csv('/Users/alexharrisson/Documents/uni/BA/output/' + i + '.tsv', dtype='unicode', sep='\t', header=0)
        for j in file.itertuples():
            #print(j[2])
            if(j[2] in hashMapOverlap):
                hashMapOverlap[j[2]] += 1
            else:
                hashMapOverlap[j[2]] = 1

    return hashMapOverlap


# build dict (hashMap) for each KB
# pos is coloum number of rsid in file
def getDict(file, pos):
    dict = {}
    for row in file.itertuples():
        rsid = str(row[pos])
        dict[rsid] = 1

    return dict


# check which ids are not the same hgvs
def checkSimilarity(dicts):
    for i in range(0, len(dicts[0])):
        # if(str(dicts[0][i][5]) == str(dicts[1][i][5])):
        print(dicts[0][i][5])


def combineBoth(dicts, kb):
    hashMapOverlap = {}
    for i in dicts:
        for j in i.items():
            if(j[0] in hashMapOverlap):
                hashMapOverlap[j[0]] += 1
            else:
                hashMapOverlap[j[0]] = 1

    file_output = open("../output/" + kb + ".tsv", "w")
    file_output.write("# \t RSID \n")

    count = 0
    for i in hashMapOverlap.items():
        count += 1;
        output = (str(count) + "\t" + str(i[0]) + "\n")
        # print(output)
        file_output.write(output)

    return hashMapOverlap

if __name__ == "__main__":
    main()
