#!/usr/bin/env python

import json
import logging
import os
from os import path
from pathlib import Path
import time
from reconstruction import reconstruct

from futils import timeit
from tqdm import tqdm

from matching import Library, Sequence, match_library

import plac


# Logging configuration
current_file = path.basename(__file__).split(".")[0]


@timeit
def match_libs(seq, libs, threshold=0.5):
    """
    Match libs with the sequence
    """
    result = []
    for lib in tqdm(libs):
        pop = {}
        pop["library"] = lib["name"]
        # See algo1_lycopene for the parameter below in the template
        # threshold = lib["score_threshold"]
        candidates = match_library(seq, Library(lib), threshold, direction53=True)
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
        result.append(pop)
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


def get_slices_libs(template):
    """
    Get slices libraries

    Args:
        template (dict): Template JSON data as a dict structure
    Returns:
        dict of slices libraries
    """

    slices_libs = {}
    for sli in template["template_slices"]:
        libs = []
        for pos in sli["template_slice"]:
            lib = template["template"]["structure"][pos - 1]["library_source"]
            libs.append(template["component_sources"][lib])
        slices_libs[sli["name"]] = libs
    return slices_libs


@timeit
def iter_all_seq(
    input_sequences,
    template_json_file,
    match_output_filename,
    reconstruction_output_filename,
    threshold=0.99,
):
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

    # Get the filenames in a list and not this freakin generator
    seq_filenames = []
    for seq in sequences:
        for filename in seq:
            seq_filenames.append(filename)

    # Loop over the sequences
    r = []
    for filename in seq_filenames:
        sq = Sequence(filename)
        json_to_output = {}
        json_to_output["target"] = sq.name

        # Logging
        logging.info(f"Target sequence: {sq.name}")

        with open(template_json_file) as json_file:
            template = json.load(json_file)

        # Get libs from template
        template["template"]
        libs = get_slices_libs(template)
        libs_to_match = libs["construct"]  # name of the fake primer

        # Match sequence
        matches = match_libs(sq, libs_to_match, threshold=threshold)
        json_to_output["matches"] = matches
        r.append(json_to_output)

    # Write output result in JSON file
    with open(match_output_filename, "w") as filename:
        json.dump(r, filename, indent=2, separators=(",", ":"))


def match(template, threshold=0.99, *targets):
    """
    Match
    """
    # Load JSON template
    with open(template) as json_file:
        template = json.load(json_file)
    r = []

    # Matching
    for target in targets:
        sq = Sequence(target)
        json_to_output = {}
        json_to_output["target"] = sq.name
        libs = get_slices_libs(template)
        libs_to_match = libs["construct"]  # name of the fake primer
        matches = match_libs(sq, libs_to_match, threshold=threshold)
        json_to_output["matches"] = matches
        r.append(json_to_output)
    s = json.dumps(r, indent=2, separators=(",", ":"))
    reconstruction_result = reconstruct(r)
    ss = json.dumps(reconstruction_result, indent=2, separators=(",", ":"))

    return ss


@plac.pos("template", "JSON construct template. Example: consruct_template.json")
@plac.pos("targets", f"Target sequence files. Example: Sanger008.seq", type=str)
@plac.opt("threshold", "Threshold", type=float)
def main(template, threshold=0.99, *targets):
    """
    cMatch command line tool
    """
    result = match(template, threshold, *targets)
    print(result)


if __name__ == "__main__":
    plac.call(main)
