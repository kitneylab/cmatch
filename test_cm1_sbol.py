import json
import logging
from os import path
from pathlib import Path
import time
from reconstruction import reconstruct

from futils import timeit
from tqdm import tqdm

from matching import Library, Sequence, match_library

TEST_DIR = "./supplementary/Testing_Algorithm_CM_1/"

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

    # Reconstruction
    reconstruction_result = reconstruct(r)
    with open(reconstruction_output_filename, "w") as filename:
        json.dump(reconstruction_result, filename, indent=2, separators=(",", ":"))


def run_test(test_params):
    """
    Run tests
    """

    targets = test_params["targets"]
    candidates = test_params["candidates"]
    nbloop = test_params["nbloop"]
    test_id = test_params["id"]
    test_name = test_params["name"]
    threshold = test_params["threshold"]

    OUTPUT_DIR = "output_results/"

    # Logging configuration
    timestr = time.strftime("%Y%m%d-%H%M%S")
    logging.basicConfig(
        format="%(asctime)s:%(levelname)s: %(message)s",
        filename=f"logs/{current_file}-{test_id}-{timestr}.log",
        encoding="utf-8",
        level=logging.DEBUG,
    )

    msg = f"Test {test_name} Target: {targets['seq_dir_names']} - Candidates: {path.basename(candidates)} - Nb runs: {nbloop}"
    logging.info(msg)

    # Iterate and match libs over the input sequences above
    for i in range(nbloop):
        msg = f"Test run: ({i+1}/{nbloop})"
        logging.info(msg)

        # Match results filename
        match_output_filename = f"{timestr}-matching-results-{current_file}-{test_id}-run-{i+1}-from-{nbloop}.json"
        match_output_filename = path.join(OUTPUT_DIR, match_output_filename)

        # Reconstruction results filename
        reconstruction_output_filename = f"{timestr}-reconstruction-{current_file}-{test_id}-run-{i+1}-from-{nbloop}.json"
        reconstruction_output_filename = output_filename = path.join(
            OUTPUT_DIR, reconstruction_output_filename
        )

        # Iterate and match libs
        iter_all_seq(
            targets,
            candidates,
            match_output_filename,
            reconstruction_output_filename,
            threshold,
        )


def test_cm1_sbol():
    """
    Test CM1 with SBOL J23101 part 
    """

    test_params = {
        "name": "CM1 SBOL",
        "id": "cm1-sbol",
        "targets": {
            "extension": ".seq",
            "sequences_path": f"{TEST_DIR}/target_sequences",
            "seq_dir_names": ["rbs_one_easy/"],
        },
        "candidates": f"{TEST_DIR}/templates/template_sbol.json",
        "nbloop": 1,
        "threshold": 0.99,
    }

    run_test(test_params)


def main():
    """
    Main
    """

    test_cm1_sbol()


if __name__ == "__main__":
    main()
