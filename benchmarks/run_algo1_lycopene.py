import json
import logging
from os import path
from pathlib import Path
import time

from futils import timeit
from tqdm import tqdm

from matching import Library, Sequence, match_library

# Logging configuration
current_file = path.basename(__file__).split(".")[0]

logging.basicConfig(
    format="%(asctime)s:%(levelname)s: %(message)s",
    filename=f"logs/{current_file}.log",
    encoding="utf-8",
    level=logging.DEBUG,
)


@timeit
def match_libs(seq, libs):
    """
    Match libs with the sequence
    """
    result = []
    for lib in tqdm(libs):
        pop = {}
        pop["library"] = lib["name"]
        threshold = lib["score_threshold"]
        candidates = match_library(seq, Library(lib), threshold, direction53=False)
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
        sequences.append(seqs) return sequences


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
def iter_all_seq(input_sequences, output_filename, templatef):
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

        with open(templatef) as json_file:
            template = json.load(json_file)

        # Get libs from template
        template["template"]
        libs = get_slices_libs(template)

        primer = sq.name.split("_")[-2]  # ex: "462006601_i08_CrtE_G06.seq"
        print(primer)
        if primer == "CrtE":
            libs_to_match = libs["revE"]
        elif primer == "CrtI":
            libs_to_match = libs["revI"]
        elif primer == "CrtB":
            libs_to_match = libs["revB"]
        else:
            raise OSError("Primer not found", sq.name)

        # Match sequence
        matches = match_libs(sq, libs_to_match)
        json_to_output["matches"] = matches
        r.append(json_to_output)

    # Write output result in JSON file
    with open(output_filename, "w") as f:
        json.dump(r, f, indent=2, separators=(",", ":"))


def run_test(test_params):
    """
    Run tests
    """

    targets = test_params["targets"]
    candidates = test_params["candidates"]
    nbloop = test_params["nbloop"]
    test_id = test_params["id"]
    test_name = test_params["name"]

    OUTPUT_DIR = "output_results/"
    OUTPUT_FILENAME = "output.json"

    msg = f"Test {test_name} Target: {targets['seq_dir_names']} - Candidates: {path.basename(candidates)} - Nb runs: {nbloop}"
    logging.info(msg)

    # Iterate and match libs over the input sequences above
    for i in range(nbloop):
        msg = f"Test run: ({i+1}/{nbloop})"
        logging.info(msg)
        timestr = time.strftime("%Y%m%d-%H%M%S")
        OUTPUT_FILENAME = f"matching-results-{current_file}-{test_id}-run-{i}-{timestr}.json"
        iter_all_seq(targets, path.join(OUTPUT_DIR, OUTPUT_FILENAME), candidates)


def test_algo1_lycopene_sanger():
    """
    Test Algo1
    All Lycopene Sanger Targets
    Lycopene inducible Template
    """

    test_params = {
        "name": "Algo1 - Lycopene - Sanger Targets - Candidate Template",
        "id": "target-lyc-sanger-template",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/repo/inducible/sequencing",
            "seq_dir_names": [
                "seq-sourcebio-21-462006601-001",
                "seq-sourcebio-25-462779601-001",
                "seq-sourcebio-28-462779601-47f67342-001",
                "seq-sourcebio-39-462743501-001",
            ],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_lycopene_sanger.json",
        "nbloop": 1,
    }

    run_test(test_params)


def test_algo1_lycopene_one():
    """
    Test Algo1
    1 Lycopene Sanger Targets
    Lycopene inducible Template
    """

    test_params = {
        "name": "Algo1 - Lycopene - Sanger Targets - Candidate Template",
        "id": "target-lyc-sanger-template",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/lyc-basic-ass-ind/",
            "seq_dir_names": ["output/1/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_lycopene.json",
        "nbloop": 1,
    }

    run_test(test_params)


def test_algo1_lycopene_12():
    """
    Test Algo1
    12 Lycopene Sanger Targets
    Lycopene inducible Template
    """

    test_params = {
        "name": "Algo1 - Lycopene - Sanger Targets - Candidate Template",
        "id": "target-lyc-sanger-template",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/lyc-basic-ass-ind/",
            "seq_dir_names": ["output/1/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_lycopene.json",
        "nbloop": 1,
    }

    run_test(test_params)


def main():
    """
    Main
    """
    #test_algo1_lycopene_sanger()
    test_algo1_lycopene_one()


if __name__ == "__main__":
    main()
