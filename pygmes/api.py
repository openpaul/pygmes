import os 
import logging
import argparse
from pygmes.exec import gmes
import shutil
import gzip
import pygmes.version  as version
from pygmes.exec import create_dir

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
        
    
def main():
    parser = argparse.ArgumentParser(description="Evaluate completeness and contamination of a MAG.")
    parser.add_argument("--input", "-i", type=str, help="path to the fasta file")
    parser.add_argument("--output", "-o", type=str, required=True, help="Path to the output folder")
    parser.add_argument("--db", "-d", type=str, required=True, help="Path to the diamond DB")
    parser.add_argument("--noclean", dest="noclean", default = True, action="store_false",required=False, help = "GeneMark-ES needs clean fasta headers and will fail if you dont proveide them. Set this flag if you don't want pygmes to clean your headers")
    parser.add_argument("--ncores", "-n", type=int, required=False, default = 1,
            help="Number of threads to use with GeneMark-ES and Diamond")
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
    
    pygmes(options.input, options.output, options.db, clean = options.noclean,
            ncores = options.ncores)
