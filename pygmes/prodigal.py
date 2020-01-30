import logging
import os
import subprocess

class prodigal:
    def __init__(self, seq, outdir, ncores):
        self.seq =seq
        self.outdir = outdir
        self.logfile = os.path.join(outdir, "prodigal.log")
        self.faa = self.run(ncores)

    def run(self, cores=1):
        logging.debug("Launching prodigal now: %s" % self.seq)
        co = os.path.join(self.outdir, "genecoord.bgk")
        faa = os.path.join(self.outdir, "prot.faa")
        lst = ["prodigal",
            "-i", self.seq,
            "-p", "meta",
           "-o", co,"-a", faa]
        try:
            # do not rerun for now if we already attempted the training once
            if not os.path.exists(faa):
                with open(self.logfile, "w") as fout:
                    subprocess.run(" ".join(lst), cwd=self.outdir, check=True, shell=True,
                                stdout = fout, stderr = fout)
            else:
                logging.info("Prodigal output already exists")
        except Exception as e:
            logging.warning("Prodigal failed on this bin")
        return(faa)

    def check_success(self):
        if os.stat(self.faa).st_size == 0:
            return False
        if not os.path.exists(self.faa):
            return False
        

        return True
