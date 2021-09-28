import concurrent.futures
import inspect
import json
import logging
from multiprocessing import Pool
from os import listdir, path
from pathlib import Path
from pprint import pprint

import psutil
from Bio import Seq, SeqIO

# Import templates
import templates.template_vio_one as template_vio_one
import templates.template_vio_ten as template_vio_ten
import templates.template_vio_20 as template_vio_20
from futils import timeit
from matching import Library, Sequence, match_library
from tqdm import tqdm
from itertools import repeat

# Get the number of CPU cores
NB_CPU_CORES = psutil.cpu_count(logical=False)


# Logging configuration
current_file = path.basename(__file__).split(".")[0]

logging.basicConfig(
    format="%(asctime)s:%(levelname)s: %(message)s",
    filename=f"logs/{current_file}.log",
    encoding="utf-8",
    level=logging.DEBUG,
)


# Algorithm Threshold between [0.0, 1.0]
MATCH_THRESHOLD = 0.99


# Define worker for concurent execution
def worker(filename):
    """
    Worker for Algo 0
    """

    # TODO change this shit
    # template =  template_vio_one
    #template = template_vio_ten
    template = template_vio_20

    sq = Sequence(filename)
    d = {}
    d["target"] = sq.name

    # Logging
    # logging.info(inspect.getframeinfo(frame).function)
    logging.info(sq.name)

    # TEMPLATE
    sub_template_libs = get_sub_templates_libs(template)
    print(sub_template_libs)
    primer_lib = "construct"
    libs_to_match = sub_template_libs[primer_lib]

    # Match sequence
    m = match_libs(sq, libs_to_match, threshold=MATCH_THRESHOLD)
    d["matches"] = m
    return d


@timeit
def match_libs(seq, libs, threshold=0.5):
    """
    Match libs with the sequence
    """
    d = []
    for lib in tqdm(libs):
        pop = {}
        pop["library"] = lib.name
        print(lib.name)
        candidates = match_library(seq, lib, threshold)
        cl = []
        for candidate in candidates:
            for c in candidate:
                cl.append(
                    {
                        "name": c.name,
                        "score": c.score,
                        "start": c.start,
                        "length": c.length,
                        "end": c.end,
                    }
                )
        pop["candidates"] = cl
        d.append(pop)
    result = d
    return result


def get_sequences(dir_dict):
    """
    Return list of sequences
    """

    SEQUENCES_EXTENSION = dir_dict["extension"]
    SEQUENCES_PATH = dir_dict["sequences_path"]
    seq_dir_names = dir_dict["seq_dir_names"]

    sequences = []
    for seq_dir in seq_dir_names:
        seqs = Path(path.join(SEQUENCES_PATH, seq_dir)).rglob(
            "*{0}".format(SEQUENCES_EXTENSION)
        )
        sequences.append(seqs)
    return sequences


def get_sub_templates_libs(template):
    """
    Get sub-templates Libraries
    """

    sublibs = {}
    for sub in template.sub_templates:
        libs = []
        for pos in sub["template_slice"]:
            lib = template.template["structure"][pos - 1]["library_source"]
            libs.append(Library(eval("template." + lib)))
        sublibs[sub["name"]] = libs
    return sublibs


@timeit
def iter_all_seq(input_sequences, output_filename, template):
    """
    Iterate over sequences

    Args:
        input_sequences (dict): Input dictionary with info about the input sequences:
        output_filename (str): Output filename

    Example:

    input_sequences = {
        'extension' = ".seq"
        'sequences_path' = "/data/Imperial/src/lyc-basic-ass-ind/"
        'seq_dir_names' = ["output"]
    }
    """

    # Get sequences to match
    sequences = get_sequences(input_sequences)
    print(sequences)

    # Get the filenames in a list and not this freakin generator
    seq_filenames = []
    for seq in sequences:
        for filename in seq:
            seq_filenames.append(filename)
    print(seq_filenames)

    # Loop over the sequences
    r = []

    print(template)
    print(seq_filenames)
    # Multiprocessing
    with Pool(NB_CPU_CORES) as p:
        # vals = p.map(worker, sorted(seq_filenames), [template for x in seq_filenames])
        vals = p.map(worker, sorted(seq_filenames))
        print(vals)
        # vals = p.imap_unordered(worker, sorted(seq_filenames))
        for d in vals:
            r.append(d)

    # Write output result in JSON file
    with open(output_filename, "w") as f:
        json.dump(r, f, indent=2, separators=(",", ":"))


def run_test(targets, candidates, nbloop):
    """
    Run tests
    """

    OUTPUT_DIR = "output_results/"
    OUTPUT_FILENAME = "output.json"

    msg = f"/// Algo0 - {targets} - {candidates} - {nbloop}"
    logging.info(msg)

    # Iterate and match libs over the input sequences above
    print(targets)
    for i in range(nbloop):
        iter_all_seq(targets, path.join(OUTPUT_DIR, OUTPUT_FILENAME), candidates)


def test_algo0_1():
    """
    Test Algo0
    """

    targets = {
        "extension": ".seq",
        "sequences_path": "/data/Imperial/src/violacein-basic-ass",
        "seq_dir_names": ["output/one/"],
    }

    candidates = template_vio_one
    nbloop = 3

    run_test(targets, candidates, nbloop)


def test_algo0_10():
    """
    Test Algo0
    """

    targets = {
        "extension": ".seq",
        "sequences_path": "/data/Imperial/src/violacein-basic-ass",
        "seq_dir_names": ["output/ten/"],
    }

    candidates = template_vio_ten
    nbloop = 3
    run_test(targets, candidates, nbloop)


def main():
    """
    Main
    """
    #test_algo0_1()
    test_algo0_10()


if __name__ == "__main__":
    main()
