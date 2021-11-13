import json
import logging
from os import path
from pathlib import Path

from futils import timeit
from tqdm import tqdm

from matching import Library, Sequence, match_library

TEST_DIR = "./supplementary/Testing_Algorithm_CM_1/"

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


@timeit
def match_libs(seq, libs, threshold=0.5):
    """
    Match libs with the sequence
    """
    result = []
    for lib in tqdm(libs):
        pop = {}
        pop["library"] = lib["name"]
        candidates = match_library(seq, Library(lib), threshold)
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
        libs_to_match = libs["construct"]  # name of the fake primer

        # Match sequence
        matches = match_libs(sq, libs_to_match, threshold=MATCH_THRESHOLD)
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
        OUTPUT_FILENAME = f"matching-results-{current_file}-{test_id}-run-{i}.json"
        iter_all_seq(targets, path.join(OUTPUT_DIR, OUTPUT_FILENAME), candidates)


def test_algo0_target_1_candidate_1():
    """
    Test Algo0
    1 Target
    1 Candidate
    """

    test_params = {
        "name": "Algo0 - 1 Target - 1 Candidate",
        "id": "target-1-candidate-1",
        "targets": {
            "extension": ".seq",
            "sequences_path": f"{TEST_DIR}/target_sequences",
            "seq_dir_names": ["target1/"],
        },
        "candidates": f"{TEST_DIR}/templates/template_vio_one.json",
        #    "nbloop": 10,
        "nbloop": 1,
    }

    run_test(test_params)


def test_algo0_target_1_candidate_1_long():
    """
    Test Algo0
    1 Target
    1 Candidate
    """

    test_params = {
        "name": "Algo0 - 1 Target - 1 Candidate",
        "id": "target-1-candidate-1",
        "targets": {
            "extension": ".seq",
            "sequences_path": f"{TEST_DIR}/target_sequences",
            "seq_dir_names": ["target1/"],
        },
        "candidates": f"{TEST_DIR}/templates/template_vio_one.json",
        #    "nbloop": 10,
        "nbloop": 1,
    }

    run_test(test_params)


def test_algo0_target_10_candidate_1():
    """
    Test Algo0
    10 Targets
    1 Candidate
    """

    test_params = {
        "name": "Algo0 - 10 Targets - 1 Candidate",
        "id": "target-10-candidate-1",
        "targets": {
            "extension": ".seq",
            "sequences_path": f"{TEST_DIR}/sequences",
            "seq_dir_names": ["./"],
            "sequences_path": f"{TEST_DIR}/target_sequences",
            "seq_dir_names": ["target10/"],
        },
        "candidates": f"{TEST_DIR}/templates/template_vio_one.json",
        #    "nbloop": 10,
        "nbloop": 1,
    }

    run_test(test_params)


def test_algo0_target_1_candidate_10():
    """
    Test Algo0
    1 Target
    10 Candidates
    """

    test_params = {
        "name": "Algo0 - 1 Target - 10 Candidates",
        "id": "target-1-candidate-10",
        "targets": {
            "extension": ".seq",
            "sequences_path": f"{TEST_DIR}/target_sequences",
            "seq_dir_names": ["target1/"],
        },
        "candidates": f"{TEST_DIR}/templates/template_vio_ten.json",
        #    "nbloop": 10,
        "nbloop": 1,
    }

    run_test(test_params)


def test_algo0_target_1_candidate_100():
    """
    Test Algo0
    1 Target
    100 Candidates
    """

    test_params = {
        "name": "Algo0 - 1 Target - 100 Candidates",
        "id": "target-1-candidate-100",
        "targets": {
            "extension": ".seq",
            "sequences_path": f"{TEST_DIR}/target_sequences",
            "seq_dir_names": ["target1/"],
        },
        "candidates": f"{TEST_DIR}/templates/template_vio_100.json",
        "nbloop": 1,
    }

    run_test(test_params)


def test_cm1_cat2min():
    """
    Test CM1
    1 Target
    1 Candidates
    """

    test_params = {
        "name": "CM1 Violacein cat2x min RBS",
        "id": "cm1-vio-cat2-min",
        "targets": {
            "extension": ".seq",
            "sequences_path": f"{TEST_DIR}/target_sequences",
            "seq_dir_names": ["targetcatx2/"],
        },
        "candidates": f"{TEST_DIR}/templates/template_vio_cat2min.json",
        "nbloop": 3,
    }

    run_test(test_params)


def main():
    """
    Main
    """
    test_cm1_cat2min()


if __name__ == "__main__":
    main()
