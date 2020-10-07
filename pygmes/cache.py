import hashlib
import logging
import os
import shutil


def compute_hash(path):
    """
    compute hash of a given file
    """
    with open(path, "rb") as f:
        bytes = f.read()  # read entire file as bytes
        return hashlib.sha256(bytes).hexdigest()


class cache:
    def __init__(self, location, fasta):
        self.location = location
        self.fasta = fasta
        self.hash = compute_hash(self.fasta)
        self.path = self._construct_path()

        return

    def _construct_path(self):
        hsh = self.hash
        return os.path.join(self.location, hsh[0:2], hsh[2:4], hsh)

    def exists(self):
        if os.path.exists(self.path):
            # should contain a hashsum file
            expected_files = [os.path.join(self.path, f) for f in ["sha256.txt", "proteins.faa", "locations.bed"]]
            for path in expected_files:
                if not os.path.exists(path):
                    return False

            # only condition to return true
            logging.debug("cache contains needed files")
            return True

        return False

    def store(self, faa, bed):
        """
        """
        if not os.path.exists(self.path):
            os.makedirs(self.path, exist_ok=True)

        # now copy files and make the sha file
        target_faa = os.path.join(self.path, "proteins.faa")
        target_bed = os.path.join(self.path, "locations.bed")
        target_sha = os.path.join(self.path, "sha256.txt")
        shutil.copy2(faa, target_faa)
        shutil.copy2(bed, target_bed)

        src_lst = [self.fasta, faa, bed]
        dst_lst = [self.fasta, target_faa, target_bed]

        # check that the copy worked and then safe the hash
        with open(target_sha, "w") as fout:
            for src, dest in zip(src_lst, dst_lst):
                h = compute_hash(src)
                h2 = compute_hash(dest)
                if h != h2:
                    logging.error("Could not store the files successfully")
                    raise OSError("Copy and source file do not have the same hash!")

                n = os.path.basename(dest)
                fout.write("{}\t{}\n".format(n, h))
        logging.debug("Wrote cache to: {}".format(self.path))
        return True

    def restore(self, target, what=None):
        """
        Restoreing function
        """
        if what is None:
            what = target[-3:]

        logging.debug("Trying to restore {} from: {}".format(what, self.path))

        if what not in ["faa", "bed"]:
            raise ValueError("what needs to be faa or bed")

        if os.path.exists(target):
            logging.warning("Wil not overwrite existsing output files")
            return

        if self.exists():
            if what == "faa":
                p = os.path.join(self.path, "proteins.faa")
            elif what == "bed":
                p = os.path.join(self.path, "locations.bed")
            shutil.copy2(p, target)
            logging.debug("Restored {}".format(p))
        else:
            return None
