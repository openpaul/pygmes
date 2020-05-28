import logging
import os
import subprocess
from pyfaidx import Fasta
from random import sample
from collections import defaultdict
from ete3 import NCBITaxa



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
        if ncores == 1:
            logging.warning("You are running Diamond with a single core. This will be slow. We recommend using 8-16 cores.")
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
        logging.debug("Inferring the lineage")
        proteinlngs = self.lineage_infer_protein(self.result)
        self.lineage  = self.vote_bin(proteinlngs)
        logging.debug("Finished the diamond step")

    def search(self, outfile, query):
        if not os.path.exists(outfile) or os.stat(outfile).st_size == 0:
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
            logging.debug("Ran diamond")
        else:
            logging.info("Diamond output already exists ")
            logging.debug("AT: %s" %  outfile)

    def sample(self, output, n=200):
        logging.debug("Sampeling %d proteins from %s" % (n, self.faa))
        try:
            faa = Fasta(self.faa)
        except ZeroDivisionError:
            logging.warning(
                "Could not read the faa file as it probably \n contains no sequence information. \n Check file: %s "
                % self.faa
            )
            return 0
        except Exception as e:
            print(e)
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
            ncbi = NCBITaxa()
            try:
                self.lineages[tax] = ncbi.get_lineage(tax)
                return self.lineages[tax]
            except ValueError:
                print(f"Not able to fetch lineage for taxid {tax}")
                return []

    def lineage_infer_protein(self, result):
        prot = {}
        for protein, taxids in result.items():
            lngs = []
            for taxid in taxids:
                l = self.inferlineage(taxid)
                if len(l) > 0:
                    lngs.append(l)

            prot[protein] = majorityvote(lngs)
        return prot

    def vote_bin(self, proteinlngs):
        lngs = [lng for prot, lng in proteinlngs.items()]
        return majorityvote(lngs)

    
class multidiamond(diamond):
    def __init__(self,proteinfiles, names, outdir, db, ncores = 1, nsample = 200):
        self.outdir = os.path.abspath(outdir)
        self.files = proteinfiles
        self.names = names
        self.samplefile = os.path.join(outdir, "samplefile.faa")
        self.outfile = os.path.join(outdir, "diamond.result")
        self.log = os.path.join(self.outdir, "diamond.log")
        self.db = db
        self.ncores = ncores
        self.lineages = {}
        if ncores == 1:
            logging.warning("You are running Diamond with a single core. This will be slow. We recommend using 8-16 cores.")
        # sample
        with open(self.samplefile, "w") as f:
            f.write("")
        for fasta, name in zip(self.files, self.names):
            self.sample(fasta, name, self.samplefile)
        
        # then run 
        self.search(self.outfile, self.samplefile)
        self.result = self.parse_results(self.outfile)
        self.lngs = self.vote_bins(self.result)

    def sample(self, fasta, name, output, n=200):
        logging.debug("Sampeling %d proteins from %s" % (n, fasta))
        try:
            faa = Fasta(fasta)
        except ZeroDivisionError:
            logging.warning(
                "Could not read the faa file as it probably \n contains no sequence information. \n Check file: %s "
                % fasta
            )
            return 0
        except Exception as e:
            print(e)
            return 0

        keys = faa.keys()
        if len(keys) > n:
            keys = sample(keys, n)
        with open(output, "a") as fout:
            for k in keys:
                fout.write(f">{name}_binseperator_{k}\n{str(faa[k])}\n")
    
    def parse_results(self, result):
        def subdict():
            return(defaultdict(list))
        r = defaultdict(subdict)
        with open(result) as f:
            for line in f:
                l = line.strip().split("\t")
                names = l[0].split("_binseperator_")
                binname = names[0]
                protein = names[1]
                r[binname][protein].append(l[5])
        return r
    
    def vote_bins(self, result):
        binnames = result.keys()
        lngs = {}
        for bin in binnames:
            protlng = self.lineage_infer_protein(result[bin])
            lngs[bin] = {}
            lngs[bin]["lng"] = self.vote_bin(protlng)
            lngs[bin]["n"] = len(protlng)

        return(lngs)
