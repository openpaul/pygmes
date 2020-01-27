import os 
import logging
import argparse
from pygmes.exec import gmes
from pygmes.diamond import multidiamond
import shutil
import gzip
from glob import glob
import pygmes.version  as version
from pygmes.exec import create_dir
from pygmes.printlngs import write_lngs

this_dir, this_filename = os.path.split(__file__)
MODELS_PATH = os.path.join(this_dir, "data", "models")




class pygmes:
    """
    Main class exposing the functionality

    Parameters:

    **fasta:** path to a fasta file

    **outdir:** path to a writable directory

    **db:** path to a diamond database with tax information

    **clean:** bool indicating if faster needs cleaning of headers

    **ncores:** number of threads to use
    """
    def __init__(self, fasta, outdir, db,  clean = True, ncores = 1):
        self.fasta = fasta
        self.outdir = outdir
        self.ncores = ncores

        if clean:
            # copy and clean file
            self.cleanfasta = self.clean_fasta(self.fasta, self.outdir)
        else:
            self.cleanfasta = self.fasta

        logging.info("Launching GeneMark-ES")
        g = gmes(self.cleanfasta, outdir, ncores)
        logging.debug("Run complete launch")
        g.run_complete(MODELS_PATH, db)
        if g.finalfaa:
            logging.debug("Copying final faa from: %s" % g.finalfaa)
            shutil.copy(g.finalfaa, os.path.join(self.outdir, "predicted_proteins.faa"))
            g.gtf2bed(g.finalgtf, os.path.join(self.outdir, "predicted_proteins.bed"))
            g.writetax()
        
    def clean_fasta(self, fastaIn, folder):
        create_dir(folder)
        name = os.path.basename(fastaIn)
        fastaOut = os.path.join(folder, name)
        mappingfile = os.path.join(folder, "mapping.csv")
        logging.debug("Cleaning fasta file")
        nms = []
        with open(fastaOut, "w") as o, open(mappingfile, "w") as mo:
            mo.write("old,new\n")
            if fastaIn.endswith(".gz"):
                openMethod = gzip.open
                gz = True
                logging.debug("reading gzipped file")
            else:
                openMethod = open
                gz = False
            # read in the fasta
            with openMethod(fastaIn) as f:
                for line in f:
                    if gz:
                        line = line.decode()
                    if line.startswith(">"):
                        line = line.strip()
                        l = line.split()
                        # get first element, usually a chromosome
                        N = l[0].strip()
                        # while the name is already taken, count up and reformat
                        n = N
                        i = 0
                        while n in nms:
                            n = "{}.{}".format(N, i)
                            i += 1
                        o.write("{}\n".format(n))
                        # finally add to know  name list
                        nms.append(n)
                        mappingline = "{},{}\n".format(line,n).replace(">", "")
                        mo.write(mappingline)
                    else:
                        o.write(line)
        return fastaOut
        
class metapygmes(pygmes):
    """
    run pygmes in metagenomic mode. This means
    we will first try the self training mode on
    all bins in the given folder
    We will then use the models from all runs
    to predict proteins in the remaining bins
    and choose the protein prediction with the largest
    number of AA. We then infer the lineage of each bin
    """
    def __init__(self, bindir, outdir, db, clean = True, ncores = 1, infertaxonomy = True):
        # find all files and 
        outdir = os.path.abspath(outdir)
        self.outdir = outdir
        bindir = os.path.abspath(bindir)
        fa = glob(os.path.join(bindir, "*.fa"))
        fna = glob(os.path.join(bindir, "*.fna"))
        fasta = glob(os.path.join(bindir, "*.fasta"))
        files = fa + fna + fasta
        # convert all files to absolute paths
        files = [os.path.abspath(f) for f in files]
        names = [os.path.basename(f) for f  in files]
        outdirs = [os.path.join(outdir, "gmes", name) for name in names]
        proteinfiles = []
        proteinnames = []
        # if needed, we clean the fasta files
        if clean:
            files = [self.clean_fasta(f, o) for f,o in zip(files, outdirs) ]

        # run self training step on all bins
        gmesm = {}
        for path, tdir, name in zip(files, outdirs, names):
            gmesm[name] = gmes(path, tdir, ncores)
            gmesm[name].selftraining()
            gmesm[name].succeed = gmesm[name].check_success()

        # this was the first round
        # now we copyall models
        modeldir = os.path.join(outdir, "1_models")
        create_dir(modeldir)
        nmodels = 0
        for name, g  in gmesm.items():
            if g.succeed:
                proteinfiles.append(g.protfaa)
                proteinnames.append(name)
                # look for the model file and if present copy into the model directory
                expectedmodel = os.path.join(g.outdir, "output","gmhmm.mod")
                if os.path.exists(expectedmodel):
                    shutil.copy(expectedmodel, os.path.join(modeldir, "{}.mod".format(name)))
                    nmodels += 1
        # now for each bin which we did not predict yet, we use the models from the dir
        if nmodels == 0:
            logging.warning("Could create models for a single bin. Thus we stop here")
            return

        # run step 2 of using the models
        for name, g in gmesm.items():
            if not g.succeed:
                g.premodel(modeldir)
                if g.bestpremodel is not False and  g.bestpremodel.check_success():
                    shutil.copy(g.bestpremodel.protfaa, g.outdir)
                    shutil.copy(g.bestpremodel.gtf, g.outdir)
                    g.succeed = True
                    proteinfiles.append(g.protfaa)
                    proteinnames.append(name)

        # diamond is faster when using more sequences
        # thus we pool all fasta together and seperate them afterwards
        diamonddir = os.path.join(outdir, "diamond")
        create_dir(diamonddir)
        logging.info("Predicting the lineage")
        dmnd = multidiamond(proteinfiles, proteinnames, diamonddir, db = db, ncores = ncores)
        logging.debug("Ran diamond and inferred lineaged")
        lineagefile = os.path.join(self.outdir, "lineages.tsv")
        write_lngs(dmnd.lngs, lineagefile)


def main():
    parser = argparse.ArgumentParser(description="Evaluate completeness and contamination of a MAG.")
    parser.add_argument("--input", "-i", type=str, help="path to the fasta file, or in metagenome mode path to bin folder")
    parser.add_argument("--output", "-o", type=str, required=True, help="Path to the output folder")
    parser.add_argument("--db", "-d", type=str, required=True, help="Path to the diamond DB")
    parser.add_argument("--noclean", dest="noclean", default = True, action="store_false",required=False, help = "GeneMark-ES needs clean fasta headers and will fail if you dont proveide them. Set this flag if you don't want pygmes to clean your headers")
    parser.add_argument("--ncores", "-n", type=int, required=False, default = 1,
            help="Number of threads to use with GeneMark-ES and Diamond")
    parser.add_argument("--meta", dest="meta", action = "store_true", default="False", help = "Run in metaegnomic mode")
    parser.add_argument(
        "--quiet", "-q", dest="quiet", action="store_true", default=False, help="Silcence most output",
    )
    parser.add_argument(
        "--debug", action="store_true", default=False, help="Debug and thus ignore safety",
    )
    parser.add_argument("-v", "--version", action="version", version=f"pygmes version {version.__version__}")
    options = parser.parse_args()

    # define logging
    logLevel = logging.INFO
    if options.quiet:
        logLevel = logging.WARNING
    elif options.debug:
        logLevel = logging.DEBUG
    logging.basicConfig(
        format="%(asctime)s %(message)s", datefmt="%m/%d/%Y %H:%M:%S: ", level=logLevel,
    )

    # check if input is readable
    if not os.path.exists(options.input):
        logging.warning("Input file does not exist: %s" % options.input)
        exit()
    logging.info("Starting pygmes")
    logging.debug("Using fasta: %s" % options.input)
    logging.debug("Using %d threads" % options.ncores)

    if not options.meta:
        pygmes(options.input, options.output, options.db, clean = options.noclean,
            ncores = options.ncores)
    else:
        metapygmes(options.input, options.output, options.db, clean = options.noclean,
            ncores = options.ncores)

