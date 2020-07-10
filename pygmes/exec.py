import logging
import os
import subprocess
import glob
import re
from pyfaidx import Fasta
from pyfaidx import FastaIndexingError
from random import sample, seed
from collections import defaultdict
from pygmes.diamond import diamond
from pygmes.printlngs import print_lngs
from ete3 import NCBITaxa
import shutil
import gzip
import urllib.request

seed(4145421)

url = "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pygmes/latest/"


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which
    from shutil import which

    return which(name) is not None


def check_dependencies(software):
    for p in software:
        if is_tool(p) is False:
            logging.error("Dependency {} is not available".format(p))
            exit(1)


def create_dir(d):
    if not os.path.isdir(d):
        try:
            os.makedirs(d)
        except OSError as e:
            logging.warning(f"Could not create dir: {d}\n{e}")


def delete_folder(d):
    if os.path.exists(d):
        if os.path.isdir(d):
            try:
                shutil.rmtree(d)
            except Exception as e:
                logging.warning("Could not delete folder: %s" % d)
                print(e)


def touch(fname, mode=0o666, dir_fd=None, **kwargs):
    flags = os.O_CREAT | os.O_APPEND
    with os.fdopen(os.open(fname, flags=flags, mode=mode, dir_fd=dir_fd)) as f:
        os.utime(
            f.fileno() if os.utime in os.supports_fd else fname, dir_fd=None if os.supports_fd else dir_fd, **kwargs
        )


class multistep_gmes:
    def __init__(self, fasta, outdir, ncores, diamonddb, models):
        self.outdir = outdir
        self.training_dir = os.path.join(self.outdir, "training")
        self.model_dir = os.path.join(self.outdir, "prediction")
        self.fasta = fasta
        self.success = False

        # call a self training session
        training = gmes(self.fasta, self.training_dir, ncores)
        training.selftraining()
        # check if training was successfulle
        if training.check_success():
            self.designate_winner(training, training=True)
        else:
            # run models to find best model
            logging.info("Using pre-trained models")
            model_run = gmes(self.fasta, self.model_dir, ncores)
            model_run.fetchinfomap()
            premodel_1 = model_run.premodel(models)
            if premodel_1 is not False:
                premodel_1.estimate_tax(diamonddb)
                print_lngs(model_run.modelinfomap[premodel_1.modelname], premodel_1.tax)
                localmodels = model_run.infer_model(premodel_1.tax)
                premodel_2 = model_run.premodel(localmodels, stage=2, existing_prediction=premodel_1)

                if premodel_2 is not False and premodel_2.check_success():
                    self.designate_winner(premodel_2)
                    # print lineage of model compared to the infered tax
                    print_lngs(model_run.modelinfomap[premodel_2.modelname], premodel_1.tax)
                else:
                    logging.error("We could not predict any proteins using pretrained models")

    def designate_winner(self, run, training=False):
        self.gmes = run
        # if training was done, copy the model
        if training:
            logging.debug("Trying to copy the model file")
            model_out = os.path.join(self.outdir, "gmhmm.mod")
            if os.path.exists(run.model):
                shutil.copy(run.model, model_out)
            else:
                logging.debug("Could not find model: {}".format(run.model))

        logging.debug("Copying final files")
        # copy the faa and bed file
        shutil.copy(run.finalfaa, os.path.join(self.outdir, "predicted_proteins.faa"))
        run.writetax()
        shutil.copy(run.bedfile, os.path.join(self.outdir, "predicted_proteins.bed"))

        self.success = True

    def cleanup(self):
        delete_folder(self.training_dir)
        delete_folder(self.model_dir)


class gmes:
    def __init__(self, fasta, outdir, ncores=1):
        self.fasta = os.path.abspath(fasta)
        self.outdir = os.path.abspath(outdir)
        self.logfile = os.path.join(self.outdir, "pygmes.log")
        self.loggtf = os.path.join(self.outdir, "pygmes_gtf.log")
        # make sure the output folder exists
        create_dir(self.outdir)
        self.ncores = ncores

        self.gtf = os.path.join(self.outdir, "genemark.gtf")
        self.protfaa = os.path.join(self.outdir, "prot_seq.faa")
        self.finalfaa = False
        self.finalgtf = False
        self.bedfile = False
        self.tax = []
        self.modelinfomap = {}
        if ncores == 1:
            logging.warning(
                "You are running GeneMark-ES with a single core. This will be slow. We recommend using 4-8 cores."
            )

    def selftraining(self):
        failpath = os.path.join(self.outdir, "tried_already")
        if os.path.exists(failpath):
            logging.info("Self-training skipped, as we did this before and it failed")
            return
        if os.path.exists(self.gtf):
            logging.info("GTF file already exists, skipping")
            self.gtf2faa()
            return

        logging.debug("Starting self-training")
        lst = [
            "gmes_petap.pl",
            "--v",
            "--fungus",
            "--ES",
            "--cores",
            str(self.ncores),
            "--min_contig",
            "5000",
            "--sequence",
            self.fasta,
        ]
        try:
            with open(self.logfile, "a") as fout:
                subprocess.run(" ".join(lst), cwd=self.outdir, check=True, shell=True, stdout=fout, stderr=fout)
        except subprocess.CalledProcessError:
            self.check_for_license_issue(self.logfile)
            touch(failpath)
            logging.info("GeneMark-ES in self-training mode has failed")
            return
        # predict and then clean
        self.model = os.path.join(self.outdir, "output", "gmhmm.mod")
        self.gtf2faa()
        self.clean_gmes_files()

    def clean_gmes_files(self):
        # clean if there are files to clean
        # this just keeps the foodprint lower
        rmfolders = ["run", "info", "data", "output/data", "output/gmhmm"]
        for folder in rmfolders:
            p = os.path.join(self.outdir, folder)
            delete_folder(p)

    def check_for_license_issue(self, logfile):
        # we do a quick search for 'icense' as this
        # string  is in every message regarding gmes licensing issues
        # if this string pops up, we need to inform the user
        with open(logfile) as fin:
            for line in fin:
                if "icense" in line:
                    logging.error(
                        "There are issues with your GeneMark-ES license. Please check that is is availiable and not expired."
                    )
                    exit(7)

    def prediction(self, model):
        self.model = model
        self.modelname = os.path.basename(model).replace(".mod", "")
        failpath = os.path.join(self.outdir, "tried_already")
        if os.path.exists(failpath):
            logging.info("Prediction skipped, as we did this before and it failed")
            return
        if os.path.exists(self.gtf):
            logging.debug("GTF file already exists, skipping")
            self.gtf2faa()
            return
        logging.debug("Starting prediction")
        lst = [
            "gmes_petap.pl",
            "--v",
            "--predict_with",
            model,
            "--cores",
            str(self.ncores),
            "--sequence",
            self.fasta,
        ]
        try:
            with open(self.logfile, "a") as fout:
                subprocess.run(" ".join(lst), cwd=self.outdir, check=True, shell=True, stdout=fout, stderr=fout)
        except subprocess.CalledProcessError:
            self.check_for_license_issue(self.logfile)
            logging.info("GeneMark-ES in prediction mode has failed")
            touch(failpath)
        # predict and then clean
        self.gtf2faa()
        self.clean_gmes_files()

    def gtf2faa(self):
        lst = ["get_sequence_from_GTF.pl", "genemark.gtf", self.fasta]
        if not os.path.exists(self.gtf):
            logging.debug("There is no GTF file")
            return
        if os.path.exists(self.protfaa):
            logging.debug("Protein file already exists, skipping")
        else:
            try:
                with open(self.loggtf, "a") as fout:
                    subprocess.run(" ".join(lst), cwd=self.outdir, check=True, shell=True, stdout=fout, stderr=fout)
            except subprocess.CalledProcessError:
                logging.warning("could not get proteins from gtf")
        # rename the proteins, to be compatibale with CAT
        # self.gtf2bed(self.gtf, self.bedfile)
        self.rename_for_CAT()

    def parse_gtf(self, gtf):
        """Given a gtf file from genemark es it extracts
        some information to create a bed file"""
        nre = re.compile(r'gene_id "([0-9]+_g)\";')

        def beddict():
            return {"chrom": None, "r": [], "strand": None}

        beds = defaultdict(beddict)
        with open(gtf) as f:
            for line in f:
                # skip comment lines
                if line.startswith("#"):
                    continue

                l = line.split("\t")
                # chrom = l[0].strip()
                # start = int(l[3])
                # stop = int(l[4])
                # strand = l[6]

                # regex match
                m = nre.findall(l[8])
                if m is not None:
                    name = m[0]
                else:
                    continue

                # save all in the dictonary
                beds[name]["chrom"] = l[0]
                beds[name]["r"].append(int(l[3]))
                beds[name]["r"].append(int(l[4]))
                beds[name]["strand"] = l[6]
        return beds

    def gtf2bed(self, gtf, outfile, rename=None, beds=None):
        """
        given a faa file and a gtf(genemark-es format)
        we will be able to create a bed file, which can be used
        with eukcc
        """
        if os.path.exists(outfile):
            logging.warning("Bedfile already exists, skipping")
            return

        # load gtf
        if beds is None:
            beds = self.parse_gtf(gtf)
        # check that keys() are contained
        for name, v in beds.items():
            if rename is not None:
                if name not in rename.keys():
                    logging.warning("Error creating bed file")
                    exit(1)

        # write to file
        with open(outfile, "w") as f:
            for name, v in beds.items():
                if rename is not None:
                    name = rename[name]
                vals = "\t".join([v["chrom"], str(min(v["r"])), str(max(v["r"])), v["strand"], name])
                f.write("{}\n".format(vals))

    def rename_for_CAT(self, faa=None, gtf=None):
        """
        renames the protein file
        to matche the format:
            >contigname_ORFNUMBER

            eg:
                >NODE_1_1
        """
        self.finalfaa = os.path.join(self.outdir, "prot_final.faa")
        self.bedfile = os.path.join(self.outdir, "proteins.bed")
        if os.path.exists(self.finalfaa) and os.path.exists(self.bedfile):
            logging.debug("Renamed faa exists, likely from previous run. Skipping this step")
            return
        if faa is None:
            faa = self.protfaa
        if gtf is None:
            gtf = self.gtf
        try:
            faa = Fasta(faa)
        except FastaIndexingError:
            logging.warning("Fastaindexing error")
            self.finalfaa = False
            return
        except Exception as e:
            logging.warning("Unhandled pyfaidx Fasta error")
            print(e)
            self.finalfaa = False
            return
        # load gtf
        beds = self.parse_gtf(gtf)
        orfcounter = defaultdict(int)
        # keep track of the renaming, so we can rename the bed
        renamed = {}
        logging.debug("Creating metadata for %s" % self.finalfaa)
        # parse and rename
        with open(self.finalfaa, "w") as fout:
            for record in faa:
                if record.name not in beds.keys():
                    logging.warning("The protein was not found in the gtf file:")
                    print("protein: %s" % record.name)
                    print("GTF file: %s" % gtf)
                    logging.warning("stopping here, this is a bug in pygmes or an issue with GeneMark-ES")
                    exit(1)
                contig = beds[record.name]["chrom"]
                orfcounter[contig] += 1
                # we use 1 as the first number, instead of the cool 0
                newprotname = "{}_{}".format(contig, orfcounter[contig])
                # keep track of the renaming, so we can rename the bed
                renamed[record.name] = newprotname
                fout.write(">{}\n{}\n".format(newprotname, record))
        # write renamed bed
        self.gtf2bed(self.gtf, self.bedfile, renamed, beds)

    def check_success(self):
        if self.finalfaa is False:
            return False
        if not os.path.exists(self.gtf):
            return False
        if not os.path.exists(self.finalfaa):
            return False

        # now more in detail
        # check if proteins are empty maybe
        with open(self.finalfaa) as fa:
            next(fa)
            for line in fa:
                if line.strip() == "":
                    return False
                break
        # if we can not open it, its of no use
        try:
            Fasta(self.finalfaa)
        except Exception:
            return False

        return True

    def estimate_tax(self, db):
        ddir = os.path.join(self.outdir, "diamond")
        create_dir(ddir)
        d = diamond(self.protfaa, ddir, db, sample=200, ncores=self.ncores)
        self.tax = d.lineage

    def premodel(self, models, stage=1, existing_prediction=None):
        logging.debug("On bin: %s" % self.fasta)
        logging.debug("Running the pre Model stage %d" % stage)
        logging.debug("Using model directory: %s", models)
        self.bestpremodel = False
        modelfiles = glob.glob(os.path.join(models, "*.mod"))
        # incoporate existing prediction if possible
        if existing_prediction is None:
            subgmes = []
        else:
            subgmes = [existing_prediction]
        logging.debug("Predicting proteins using {} models".format(len(modelfiles)))
        for model in modelfiles:
            logging.debug("Using model %s" % os.path.basename(model))
            name = os.path.basename(model)
            odir = os.path.join(self.outdir, "{}_premodels".format(stage), name)
            g = gmes(self.fasta, odir, ncores=self.ncores)
            g.prediction(model)
            if g.check_success():
                subgmes.append(g)
                if stage == 1:
                    logging.debug("Stopping stage 1 prediction, as we have one proteome to predict the lineage")
                    break

        if len(subgmes) == 0:
            logging.warning("Could not predict any proteins in this file")
            return False
        else:
            aminoacidcount = []
            for g in subgmes:
                i = 0
                try:
                    fa = Fasta(g.protfaa)
                    for seq in fa:
                        i += len(seq)
                except FastaIndexingError:
                    logging.warning("Could not read fasta")
                except Exception as e:
                    logging.debug("Unhandled pyfaidx Fasta error")
                    print(e)

                aminoacidcount.append(i)
            if len(aminoacidcount) == 0:
                logging.error("Could not determine best model")
                return False
            # set the best model as the model leading to the most amino acids
            idx = aminoacidcount.index(max(aminoacidcount))
            logging.info("Best model set as: %s" % os.path.basename(subgmes[idx].model))
            return subgmes[idx]

    def fetchinfomap(self):
        """
        function to make sure the information of all models
        is known to the class
        """
        ncbi = NCBITaxa()
        if len(self.modelinfomap) == 0:
            info = self.fetch_info("{}info.csv".format(url))
            logging.debug("Fetching models from {}".format(url))
            for line in info.split("\n"):
                l = line.split(",")
                # fetch lineage from ete3 for each model
                # time consuming but important to adapt to changes
                # in NCBI taxonomy
                if len(l) > 1:
                    self.modelinfomap[l[0]] = ncbi.get_lineage(l[1])

    def infer_model(self, tax, n=3):
        """
        given we infered a lineage or we know a lineage
        we can try to fetch a model from the number of
        precomputed models that already exists
        for this we choose the model that shares the most
        taxonomic element with the predicted lineage
        If multiple modles have similar fit, we just again chose the best one
        """
        self.fetchinfomap()
        logging.debug("Inferring model")

        candidates = self.score_models(self.modelinfomap, tax, at_least=n)

        if len(candidates) > n:
            candidates = sample(candidates, n)
            logging.debug("Reduced models to {}".format(n))
        # for each candidate, try to download the model into a file
        modeldir = os.path.join(self.outdir, "models")
        delete_folder(modeldir)
        create_dir(modeldir)
        for model in candidates:
            self.fetch_model(modeldir, url, model)

        return modeldir

    def fetch_model(self, folder, url, name):
        url = "{}/models/{}.mod.gz".format(url, name)
        modelfile = os.path.join(folder, "{}.mod".format(name))
        response = urllib.request.urlopen(url)
        data = response.read()  # a `bytes` object
        content = gzip.decompress(data)
        with open(modelfile, "w") as mod:
            mod.writelines(content.decode())

    def fetch_info(self, url, i=5):
        logging.debug("opening url {}".format(url))
        try:
            response = urllib.request.urlopen(url)
            data = response.read()  # a `bytes` object
            infocsv = data.decode("utf-8")
        except urllib.error.URLError:
            import time

            if i > 0:
                i = i - 1
                time.sleep(5)
                infocsv = self.fetch_info(url, i)
            else:
                logging.error("Could not fetch model file")
                exit(1)

        return infocsv

    def score_models(self, infomap, lng, at_least=3):
        logging.debug("scoring all models")
        scores = []
        candidates = []

        for model, mlng in infomap.items():
            score = len(set(lng) & set(mlng))
            scores.append((model, score))
        # sort scored models
        scores.sort(key=lambda x: x[1], reverse=True)
        if len(scores) == 0:
            logging.error("No models were obtained, so none were scored")
            exit(1)

        maxscore = scores[0][1]
        # get all models with the highest score,
        # and then fill up till at_least
        for x in scores:
            if x[1] == maxscore:
                candidates.append(x[0])
                logging.debug("Choose model {} with score {}".format(x[0], x[1]))
            elif len(candidates) < at_least and x[1] != maxscore:
                candidates.append(x[0])
                logging.debug("Choose model {} with score {}".format(x[0], x[1]))
        return candidates

    def writetax(self):
        """
        write infered taxonomy in a machine and human readble format
        """
        logging.info("Translating lineage")
        ncbi = NCBITaxa()
        taxf = os.path.join(self.outdir, "lineage.txt")
        with open(taxf, "w") as fout:
            # get the information
            lng = self.tax
            nms = ncbi.get_taxid_translator(lng)
            ranks = ncbi.get_rank(lng)
            # first line is taxids in machine readable
            s = "-".join([str(i) for i in lng])
            fout.write("#taxidlineage: {}\n".format(s))
            fout.write("taxid\tncbi_rank\tncbi_name\n")
            for taxid in lng:
                if taxid in nms.keys():
                    name = nms[taxid]
                else:
                    name = "unnamed"
                fout.write(f"{taxid}\t{ranks[taxid]}\t{name}\n")
        logging.info("Wrote lineage to %s" % taxf)
