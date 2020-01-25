import logging
import os
import subprocess
from pyfaidx import Fasta
from random import sample
from collections import defaultdict
from ete3 import NCBITaxa

ncbi = NCBITaxa()


def majorityvote(lngs, fraction=0.6):
    if fraction <= 0.5:
        logging.warning("fraction must be larger than 0.5")
    if len(lngs) == 0:
        return []
    ml = max([len(v) for v in lngs])
    i = 0
    n = len(lngs)

    lng = []
    while i < ml:
        choices = defaultdict(int)
        for l in lngs:
            if len(l) <= i:
                continue
            choices[l[i]] += 1
        # decide if this is majority material
        if len(choices) > 0:
            best = max(choices, key=lambda key: choices[key])
            if choices[best] / n >= fraction:
                lng.append(best)
            else:
                break
            i += 1
        else:
            break
    return lng


class diamond:
    def __init__(self, faa, outdir, db, ncores=1, sample=100):
        self.faa = faa
        self.outdir = outdir
        self.db = db
        self.ncores = ncores
        self.outfile = os.path.join(self.outdir, "diamond.results.tsv")
        self.log = os.path.join(self.outdir, "diamond.log")
        self.lineages = {}
        # sample n proteins
        logging.info("Subsampeling %d proteins" % sample)
        self.samplefile = os.path.join(self.outdir, "diamond.query.faa")
        self.sample(self.samplefile, sample)

        # runa search
        logging.info("Running diamond blastp")
        self.search(self.outfile, self.samplefile)
        logging.debug("Parsing diamond output")
        self.parse_results(self.outfile)
        # infer lineages
        logging.debug("Infering the lineage")
        self.lineage_infer_protein()
        self.vote_bin()
        logging.debug("Finsihed the diamond step")

    def search(self, outfile, query):
        if not os.path.exists(outfile):
            logging.info("Running diamond now")
            lst = [
                "diamond",
                "blastp",
                "--db",
                self.db,
                "-q",
                query,
                "-p",
                str(self.ncores),
                "--evalue",
                str(1e-20),
                "--max-target-seqs",
                "3",
                "--outfmt",
                "6",
                "qseqid",
                "sseqid",
                "pident",
                "evalue",
                "bitscore",
                "staxids",
                "-o",
                outfile,
            ]
            with open(self.log , "w") as fout:
                subprocess.run(lst, stderr = fout, stdout = fout)
        else:
            logging.info("Diamond output already exists ")
            logging.debug("AT: %s" %  outfile)

    def sample(self, output, n=200):
        try:
            faa = Fasta(self.faa)
        except ZeroDivisionError:
            logging.warning(
                "Could not read the faa file as it probably \n contains no sequence information. \n Check file: %s "
                % self.faa
            )
            return 0

        keys = faa.keys()
        if len(keys) > n:
            keys = sample(keys, n)
        with open(output, "a") as fout:
            for k in keys:
                fout.write(f">{k}\n{str(faa[k])}\n")

    def parse_results(self, result):
        r = defaultdict(list)
        with open(result) as f:
            for line in f:
                l = line.strip().split("\t")
                r[l[0]].append(l[5])
        self.result = r

    def inferlineage(self, tax):
        if tax in self.lineages.keys():
            return self.lineages[tax]
        else:
            try:
                self.lineages[tax] = ncbi.get_lineage(tax)
                return self.lineages[tax]
            except ValueError:
                print(f"Not able to fetch lineage for taxid {tax}")
                return []

    def lineage_infer_protein(self):
        prot = {}
        for protein, taxids in self.result.items():
            lngs = []
            for taxid in taxids:
                l = self.inferlineage(taxid)
                if len(l) > 0:
                    lngs.append(l)

            prot[protein] = majorityvote(lngs)
        self.proteinlngs = prot

    def vote_bin(self):
        lngs = [lng for prot, lng in self.proteinlngs.items()]
        self.lineage = majorityvote(lngs)
