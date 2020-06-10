import os 
import logging
import argparse
from pygmes.exec import gmes
from pygmes.diamond import multidiamond
import shutil
import gzip
from glob import glob
from pyfaidx import Fasta
import pygmes.version  as version
from pygmes.exec import create_dir
from pygmes.printlngs import write_lngs
from pygmes.prodigal import prodigal

this_dir, this_filename = os.path.split(__file__)
MODELS_PATH = os.path.join(this_dir, "data", "models")



class bin:
    def __init__(self, path, outdir):
        self.fasta = os.path.abspath(path)
        self.name = os.path.basename(path)
        self.outdir = os.path.join(os.path.abspath(outdir), self.name)
        create_dir(self.outdir)
        self.hybridfaa = None
    
    def get_best_faa(self):
        if self.kingdom is not None and self.kingdom in ["bacteria", "archaea"]:
            if self.prodigal.check_success():
                return(self.prodigal.faa, self.prodigal.bed, self.name, "prodigal")
        elif  hasattr(self, 'gmes') and self.gmes.check_success():
            # check if we made a hybrid
            if self.hybridfaa is not None and os.path.exists(self.hybridfaa):
                return(self.hybridfaa, self.hybridbed, self.name, "hybrid")
            else:
                return(self.gmes.finalfaa, self.gmes.bedfile, self.name, "GeneMark-ES")
        elif hasattr(self, "gmes") and not self.gmes.check_success():
            # this bin is euakryotic maybe, but gmes has failed. 
            # so we will return the prodigal peptides instead, eeventhough we 
            # know this might be of low quality
            if self.prodigal.check_success():
                return(self.prodigal.faa, self.prodigal.bed, self.name, "prodigal")
        # default to none, as we dont have proteins it seems
        logging.debug("No final faa for bin: %s" % self.name)
        return(None, None, self.name, None)

    def gmes_training(self, ncores = 1):
        outdir = os.path.join(self.outdir, "gmes_training")
        self.gmes = gmes(self.fasta, outdir, ncores)
        self.gmes.selftraining()
    
    def run_prodigal(self, ncores = 1, outdir=None):
        if outdir is None:
            outdir = os.path.join(self.outdir, "prodigal")
            create_dir(outdir)
        self.prodigal = prodigal(self.fasta,  outdir, ncores)

    def make_hybrid_faa(self, gmesfirst = True):
        logging.debug("Making a hybrid of bin %s" % self.name)
        try:
            if not self.gmes.check_success() or not self.prodigal.check_success():
                return
        except AttributeError:
            return
        outdir = os.path.join(self.outdir, "hybrid")
        create_dir(outdir)
        self.hybridfaa = os.path.join(outdir, "gmes_prodigal_merged.faa")
        self.hybridbed = os.path.join(outdir, "gmes_prodigal_merged.bed")

        def sane_faa(faa):
            if os.stat(faa).st_size == 0:
                return False
            try:
                fa = Fasta(faa)
                if len(fa.keys()) > 0:
                    return fa
                return False
            except Exception as e:
                print(e)
                logging.debug("Fasta has no entries: %s" % faa)
                return False

        def chromname(s):
            if s.startswith(">"):
                s = s[1:]
            return s.strip().split()[0].rsplit("_", 1)[0]

        def getchroms(fa):
            contigs = set()
            for seq in fa:
                contigs.add(chromname(seq.name))
            return contigs

        if gmesfirst:
            faa1 = self.gmes.finalfaa
            faa2 = self.prodigal.faa
            bed1 = self.gmes.bedfile
            bed2 = self.prodigal.bed
        else:
            faa2 = self.gmes.finalfaa
            faa1 = self.prodigal.faa
            bed2 = self.gmes.bedfile
            bed1 = self.prodigal.bed
        # load and check for valid fastas
        fa1  = sane_faa(faa1)
        fa2  = sane_faa(faa2)
        if fa1 is False or fa2 is False:
            return faa1

        # find contigs uniqly annotated in faa2
        contigs1 = getchroms(fa1)
        contigs2 = getchroms(fa2)
        leftover = contigs2 - contigs1
        if len(leftover) > 0:
            logging.debug("We found possible bacterial proteins in this proteome")
            # write prot from fa1
            with open(self.hybridfaa, "w") as fout:
                for seq in fa1:
                    fout.write(f">{seq.name}\n{seq}\n")
                # add new from fa2
                for seq in fa2:
                    # write to file
                    if chromname(seq.name) in leftover:
                        fout.write(f">{seq.name}\n{seq}\n")

            # make a merged bedfile
            shutil.copy(bed1, self.hybridbed)
            with open(self.hybridbed, "a") as fout, open(bed2) as fin:
                for line in fin:
                    if line.split("\t")[0] in leftover:
                        fout.write(line)

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
            g.writetax()
        if g.bedfile:
            shutil.copy(g.bedfile, os.path.join(self.outdir, "predicted_proteins.bed"))
        else:
            logging.debug("Could not find bed file")

    def clean_fasta(self, fastaIn, folder):
        create_dir(folder)
        name = os.path.basename(fastaIn)
        fastaOut = os.path.join(folder, name)
        if os.path.abspath(fastaOut) == os.path.abspath(fastaIn):
            logging.error("Name collision, please do not use the same folder for input and output")
            exit(1)
            
        if os.path.exists(fastaOut):
            logging.warning("Clean fasta file %s already exists, this could be from a previous run or an file name issue!" % name)
            return fastaOut
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
    def __init__(self, bindir, outdir, db, clean = True, ncores = 1, infertaxonomy = True, fill_bac_gaps = True):
        # find all files and 
        outdir = os.path.abspath(outdir)
        self.outdir = outdir
        bindir = os.path.abspath(bindir)
        fa = glob(os.path.join(bindir, "*.fa"))
        fna = glob(os.path.join(bindir, "*.fna"))
        fasta = glob(os.path.join(bindir, "*.fasta"))
        files = fa + fna + fasta
        prodigaldir = os.path.join(self.outdir, "prodigal")
        # convert all files to absolute paths
        files = [os.path.abspath(f) for f in files]
        names = [os.path.basename(f) for f  in files]
        if len(names) != len(set(names)):
            logging.warning("Bin files need to have unique names")
            exit(1)

        #outdirs = [os.path.join(outdir, name) for name in names]
        #proteinfiles = []
        #proteinnames = []
        # if needed, we clean the fasta files
    
        if clean:
            logging.info("Cleaning input fastas")
            cleanfastadir = os.path.join(outdir, "fasta_clean")
            files = [self.clean_fasta(f, cleanfastadir) for f in files ]


        # bin list to keep all the bins and handle all the operations
        binlst = []
        bindirs = os.path.join(outdir, "bins")
        for path in files:
            binlst.append(bin(path, bindirs))

        # run prodigal
        logging.info("Running prodigal on all bins")
        for b in binlst:
            b.run_prodigal(ncores = ncores)
        
        # now we can already get a first lineage estimation
        # diamond is faster when using more sequences
        # thus we pool all fasta together and seperate them afterwards
        diamonddir = os.path.join(outdir, "diamond", "step_1")
        create_dir(diamonddir)
        logging.info("Predicting the lineage")
        proteinfiles = [b.prodigal.faa for b in binlst if b.prodigal.check_success()]
        proteinnames = [b.name for b in binlst if b.prodigal.check_success()]
        dmnd_1 = multidiamond(proteinfiles, proteinnames, diamonddir, db = db, ncores = ncores)
        logging.debug("Ran diamond and inferred lineages")
        # assign a taxonomic kingdom based on the first lineage estimation
        anyeuks = False
        for b in binlst:
            b.first_lng_estimation = None
            b.kingdom = None
            if b.name in dmnd_1.lngs.keys():
                # as no lng was infered for this bin, we could try prodigal
                b.first_lng_estimation = dmnd_1.lngs[b.name]
                if 2 in b.first_lng_estimation['lng']:
                    b.kingdom = "bacteria"
                elif 2759 in b.first_lng_estimation['lng']:
                    b.kingdom = "eukaryote"
                    anyeuks = True
                elif 2157 in b.first_lng_estimation['lng']:
                    b.kingdom = "archaea"
                else:
                    anyeuks = True

        if anyeuks == False:
            logging.info("All bins are prokaryotes, we can skip the GeneMark-ES steps")
        else:
            # for all bins that are euakryotic or could not be assigned a lineage,
            # we try GeneMark-ES in a two step mode
            modeldir = os.path.join(outdir, "gmes_models")
            create_dir(modeldir)
            nmodels = 0
            logging.info("Running GeneMark-ES in self training")
            for b in binlst:
                if b.kingdom is None or b.kingdom == "eukaryote":
                    # run self training
                    b.gmes_training(ncores = ncores)
                    expectedmodel = os.path.join(b.gmes.outdir, "output","gmhmm.mod")
                    if os.path.exists(expectedmodel):
                        shutil.copy(expectedmodel, os.path.join(modeldir, "{}.mod".format(b.name)))
                        nmodels += 1
            # check if any bins were not predicted, if so we can use the models
            # from other bins to get a better estimate
            # if thats not possible, we could still run pygmes in non metagenomic 
            # on each bin, but that should be decied by the user
            if nmodels == 0:
                logging.debug("No models were successfully trained")
            else:
                for b in binlst:
                    if b.kingdom is None or b.kingdom == "eukaryote":
                        if b.gmes.check_success() is False:
                            b.gmes.premodel(modeldir)
                            # if successfull, overwrite the gmes, with the successfull gmes
                            if b.gmes.bestpremodel is not False and b.gmes.bestpremodel.check_success():
                                b.gmes = b.gmes.bestpremodel
            # now we have proteins predicted for all
            # we can now give each bin the chance to merge prodigal and Gmes predictions
            for b in binlst:
                b.make_hybrid_faa()

            # now we update the lineages using the new proteins
            # and then we can create a final set of protein files
            diamonddir = os.path.join(outdir, "diamond", "step_2")
            create_dir(diamonddir)
            proteinfiles = []
            proteinnames = []
            for b in binlst:
                path, bedpath, name, software = b.get_best_faa()
                if path is not None and b.kingdom not in ["bacteria", "archaea"]:
                    proteinfiles.append(path)
                    proteinnames.append(name)
            if len(proteinfiles) > 0:
                logging.info("Predicting the lineage using the results from GeneMark-ES")
                dmnd_2 = multidiamond(proteinfiles, proteinnames, diamonddir, db = db, ncores = ncores)
                for b in binlst:
                    if b.name in dmnd_2.lngs.keys():
                        # as no lng was infered for this bin, we could try prodigal
                        b.first_lng_estimation = dmnd_2.lngs[b.name]
                        if 2 in b.first_lng_estimation['lng']:
                            b.kingdom = "bacteria"
                        elif 2759 in b.first_lng_estimation['lng']:
                            b.kingdom = "eukaryote"
                        elif 2157 in b.first_lng_estimation['lng']:
                            b.kingdom = "archaea"
            else:
                logging.info("No changes after applying GeneMark-ES")

        # now we can make a final FAA folder:
        finaloutdir = os.path.join(self.outdir, "predicted_proteomes")
        finalbeddir = os.path.join(finaloutdir, "bed")
        create_dir(finaloutdir)
        create_dir(finalbeddir)
        lngs = {}
        metadataf = os.path.join(self.outdir, "metadata.tsv")
        metadata = {}
        finalfaas = {}
        for b in binlst:
            t =  os.path.join(finaloutdir, "{}.faa".format(b.name))
            bt = os.path.join(finalbeddir, "{}.bed".format(b.name))
            path, bedpath, name, software = b.get_best_faa()
            b.software = software
            metadata[b.name] = {"path": path, 
                             "software": software, 
                             "nprot": None,
                             "lng": [],
                             "name": b.name}
            if path is not None:
                shutil.copy(path, t)
                shutil.copy(bedpath, bt)
                metadata[b.name]['path'] = t
                finalfaas[b.name] = {"faa": t, "fasta": b.fasta}
                try:
                    fa = Fasta(t)
                    metadata[b.name]['nprot'] = len(fa.keys())
                except Exception as e:
                    print(e)
                    metadata[b.name]['nprot'] = 0
            if b.first_lng_estimation is not None:
                lngs[b.name] = {}
                lngs[b.name]['lng'] = b.first_lng_estimation['lng']
                lngs[b.name]['n'] = b.first_lng_estimation['n']
                metadata[b.name]['lng'] = "-".join([str(x) for x in b.first_lng_estimation['lng']])
        logging.debug("Copied files, now writing lineages")
        lngfile = os.path.join(outdir, "lineages.tsv")
        write_lngs(lngs, lngfile)
        # write metadata to disk
        logging.debug("Writing metadata")
        with open(metadataf, "w") as fout:
            keys = ['name', "path", "software", "nprot", "lng"] 
            fout.write("\t".join(keys))
            fout.write("\n")
            for k,v in metadata.items():
                l = []
                for key in keys:
                    l.append(str(v[key]))
                fout.write("\t".join(l))
                fout.write("\n")
        
        # make a single file for CAT, including a prefix so CAT will not get confused
        def single_fasta(fastas, names, output, sep = "_"):
            if len(fastas) != len(names):
                logging.warning("Number of Fastas does not match names")
                exit(1)
            nseqs = 0
            with open (output, "w") as fout:
                for fasta, name in zip(fastas, names):
                    with open(fasta) as fin:
                        for line in fin:
                            if line.startswith(">"):
                                line = line.strip().split()[0][1:]
                                line = ">{}{}{}\n".format(name, sep, line)
                                nseqs += 1
                            fout.write(line)
            if nseqs == 0:
                logging.warning("No sequence in aggregate")
                exit(1)
        
        # make massive protein file:
        catdir = os.path.join(outdir, "CAT")
        create_dir(catdir)
        catfaa = os.path.join(catdir, "cat.faa")
        catfna = os.path.join(catdir, "cat.fna")
        names = list(finalfaas.keys())
        names.sort()
        faas = [finalfaas[name]['faa'] for name in names]
        fnas = [finalfaas[name]['fasta'] for name in names]
        single_fasta(fnas, names, catfna)
        single_fasta(faas, names, catfaa)




        logging.info("Successfully ran pygmes --meta")


def main():
    parser = argparse.ArgumentParser(description="Evaluate completeness and contamination of a MAG.")
    parser.add_argument("--input", "-i", type=str, help="path to the fasta file, or in metagenome mode path to bin folder")
    parser.add_argument("--output", "-o", type=str, required=True, help="Path to the output folder")
    parser.add_argument("--db", "-d", type=str, required=True, help="Path to the diamond DB")
    parser.add_argument("--noclean", dest="noclean", default = True, action="store_false",required=False, help = "GeneMark-ES needs clean fasta headers and will fail if you dont proveide them. Set this flag if you don't want pygmes to clean your headers")
    parser.add_argument("--ncores", "-n", type=int, required=False, default = 1,
            help="Number of threads to use with GeneMark-ES and Diamond")
    parser.add_argument("--meta", dest="meta", action = "store_true", default=False, help = "Run in metaegnomic mode")
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
        exit(1)
    logging.info("Starting pygmes")
    logging.debug("Using fasta: %s" % options.input)
    logging.debug("Using %d threads" % options.ncores)

    if not options.meta:
        pygmes(options.input, options.output, options.db, clean = options.noclean,
            ncores = options.ncores)
    else:
        metapygmes(options.input, options.output, options.db, clean = options.noclean,
            ncores = options.ncores)

