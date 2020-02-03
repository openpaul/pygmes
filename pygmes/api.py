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
from pygmes.prodigal import prodigal

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
        if os.path.exists(fastaOut):
            logging.debug("Clean fasta %s exists" % name)
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
        # convert all files to absolute paths
        files = [os.path.abspath(f) for f in files]
        names = [os.path.basename(f) for f  in files]
        outdirs = [os.path.join(outdir, "gmes", name) for name in names]
        proteinfiles = []
        proteinnames = []
        # if needed, we clean the fasta files
        if clean:
            files = [self.clean_fasta(f, o) for f,o in zip(files, outdirs) ]

        # SELF-TRAINING
        gmesm = {}
        for path, tdir, name in zip(files, outdirs, names):
            gmesm[name] = gmes(path, tdir, ncores)
            gmesm[name].selftraining()
            gmesm[name].succeed = gmesm[name].check_success()

        # PREDICTION WITH LEARNED MODELS

        # Copy models into one location, so its easy to reuse
        modeldir = os.path.join(outdir, "1_models")
        create_dir(modeldir)
        nmodels = 0
        for name, g  in gmesm.items():
            if g.succeed:
                logging.debug("For bin %s we set finalfaa: %s" %(name, g.finalfaa))
                proteinfiles.append(g.finalfaa)
                proteinnames.append(name)
                # look for the model file and if present copy into the model directory
                expectedmodel = os.path.join(g.outdir, "output","gmhmm.mod")
                if os.path.exists(expectedmodel):
                    shutil.copy(expectedmodel, os.path.join(modeldir, "{}.mod".format(name)))
                    nmodels += 1

        # now for each bin which we did not predict yet, we use the models from the dir
        if nmodels == 0:
            logging.warning("Could create models for a single bin. Thus we stop here")
            exit(0)

        # run step 2 of using the models
        for name, g in gmesm.items():
            if not g.succeed:
                g.premodel(modeldir)
                if g.bestpremodel is not False and  g.bestpremodel.check_success():
                    g.finalfaa = g.bestpremodel.finalfaa
                    logging.debug("For bin %s we set finalfaa: %s" %(name, g.finalfaa))
                    shutil.copy(g.bestpremodel.finalfaa, g.outdir)
                    shutil.copy(g.bestpremodel.gtf, g.outdir)
                    g.succeed = True
                    proteinfiles.append(g.finalfaa)
                    proteinnames.append(name)

        # diamond is faster when using more sequences
        # thus we pool all fasta together and seperate them afterwards
        diamonddir = os.path.join(outdir, "diamond", "step_1")
        create_dir(diamonddir)
        logging.info("Predicting the lineage")
        print(proteinfiles)
        dmnd = multidiamond(proteinfiles, proteinnames, diamonddir, db = db, ncores = ncores)
        logging.debug("Ran diamond and inferred lineages")
        #lineagefile = os.path.join(self.outdir, "lineages.tsv")
        #write_lngs(dmnd.lngs, lineagefile)
        # now we know which mags are eukayotic
        # we could now run prodigal for all non eukayotic bins to see
        # if we can improve the proteinpreidction
        prodigaldir = os.path.join(self.outdir, "prodigal")
        create_dir(prodigaldir)
        logging.info("Trying prodigal for failed and non eukaryotic bins")
        for name, g in gmesm.items():
            g.try_prodigal = False
            if name not in dmnd.lngs.keys():
                # as no lng was infered for this bin, we could try prodigal
                g.try_prodigal = True
            elif fill_bac_gaps:
                # even for eukaryotes and everyone else, we will try to annotate 
                # prokaryotic genes, to fill in contigs with missing data
                g.try_prodigal = True
            else:
                # if bin was not classified as eukaryotic we try prodigal
                if 2759 not in dmnd.lngs[name]['lng']:
                    g.try_prodigal = True
            if g.try_prodigal:
               g.prodigaldir = os.path.join(prodigaldir, name)
               create_dir(g.prodigaldir)
               g.prodigal = prodigal(g.fasta,  g.prodigaldir, ncores)
            
        if fill_bac_gaps:
            logging.info("Will try to see if prodigal found proteins on contigs that GeneMark-ES missed")
            def zip_faa(faa1, faa2, outfile):
                def sane_faa(faa):
                    if os.stat(faa).st_size == 0:
                        return False
                    try:
                        fa = Fasta(faa)
                        if len(fa) > 0:
                            return fa
                        return False
                    except:
                        return False
                def getchroms(fa):
                    contigs = set()
                    for seq in fa:
                        cn = seq.name.split()[0]
                        l = cn.split("_")
                        chrom = "".join(l[0:-1])
                        contigs.add(chrom)
                    return contigs
                # load and check for valid fastas
                fa1  = sane_faa(faa1)
                fa2  = sane_faa(faa2)
                if fa is False or fa2 is False:
                    return faa1
                # find contigs uniqly annotated in faa2
                contigs1 = getchroms(fa1)
                contigs2 = getchroms(fa2)
                leftover = contigs2 - contigs1
                if len(leftover) > 0:
                    logging.debug("We found possible bacterial proteins in this proteome")
                    # write prot from fa1
                    with open(outfile, "w") as fout:
                        for seq in fa1:
                            fout.write(f">{seq.name}\n{seq}\n")
                        # add new from fa2
                        for seq in fa2:
                            cn = seq.name.split()[0]
                            l = cn.split("_")
                            chrom = "".join(l[0:-1])
                            # write to file
                            if chrom in leftover:
                                fout.write(f">{seq.name}\n{seq}\n")
                    return outfile
                else:
                    return faa1

            # some eukaryotic bins might contain single or multiple
            # bac contigs, we by comparing prodigal to GeneMark-ES proteome
            # predictins
            for name, g in gmesm.items():
                # if Gmes and prodigal worked
                if g.succeed and g.prodigal.check_success():
                    outfile = os.path.join(g.outdir,"hybrid_gmes_prodigla.faa")
                    faa_hybrid = zip_faa(g.finalfaa, g.prodigal.faa, outfile)
                    # todog.finalfaa
                    if faa_hybrid != outfile:
                        g.finalfaa = faa_hybrid

        # summarise the final result
        finaloutdir = os.path.join(self.outdir, "hybridpredicted")
        create_dir(finaloutdir)
        proteinfiles  = []
        proteinnames = []
        for name, g in gmesm.items():
            logging.debug("Copy final files: %s" % name)
            t =  os.path.join(finaloutdir, "{}.faa".format(name))
            if not g.try_prodigal:
                # all likely eukryotic bins
                print(f"cp {g.finalfaa} {t}")
                if g.succeed or (g.bestpremodel is not False and  g.bestpremodel.check_success()):
                    shutil.copy2(g.finalfaa, t)
                    proteinfiles.append(t)
                    proteinnames .append(name)
                else:
                    logging.info("Could not predict proteins for bin %s" % name)
            else:
                # all possible bacterial bins
                if os.path.exists(g.prodigal.faa) and g.check_success():
                    shutil.copy2(g.prodigal.faa, t)
                    proteinfiles.append(t)
                    proteinnames .append(name)
                else:
                    logging.info("Could not predict proteins for bin %s" % name)

        diamonddir = os.path.join(outdir, "diamond", "step_2")
        create_dir(diamonddir)
        logging.info("Predicting the lineage")
        dmnd = multidiamond(proteinfiles, proteinnames, diamonddir, db = db, ncores = ncores)
        logging.debug("Ran diamond and inferred lineages")
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

