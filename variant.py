import re, sys


class Variant:

    def __init__(self):
        """Initialize a new variant."""
        self.hgvs = None  # full hgvs notation
        self.valid = False  # check if hgvs notation is valid
        self.hgvs_error = None  # check if mutation on gene exists

        self.ac = None  # accession number from hgvs found in KB
        self.posedit = None  # mutation from hgvs found in KB

        self.entrezid = None  # Entrez Gene Id/Locus ID

        self.assembly = None  # assebly (eg. 37 or 38)

        self.refseq = list()  # RefSeq ID
        self.ens = None  # Ensemble ID

        self.rsid = -1  # dbSNP ID if found (in KB)

        self.mydb = None  # Used for mapping ens->refseq
        self.mycursor = None  # Used for mapping ens->refseq

    def printVariant(self):
        """Print Variant information."""
        print(
            'found HGVS: ' + str(self.hgvs) + '\n'
            'is valid HGVS: ' + str(self.valid) + '\n'
            'used assembly: ' + str(self.assembly) + '\n'
            'Entrez ID: ' + str(self.entrezid) + '\n'
            'RefSeq ID: ' + str(self.refseq) + '\n'
            'Ensemble ID: ' + str(self.ens) + '\n'
            'RSID: ' + str(self.rsid) + '\n')

    def proteinOneLetter(self):
        """Convert 1 letter amino acids to 3 letter."""
        aminoAcids = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
                      'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
                      'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W',
                      'Ala': 'A', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'}

        d = {v: k for k, v in aminoAcids.items()}  # invert dictionary

        try:
            match = re.match(r"([a-z]+)([0-9]+)([a-z]+)", self.posedit, re.I)
            items = list()
            if match:
                items = match.groups()
                self.posedit = "p." + d[items[0]] + items[1] + d[items[2]]
        except KeyboardInterrupt:
            print("Exiting.. \n")
            sys.exit(0)
        except Exception:
            self.posedit = ""
