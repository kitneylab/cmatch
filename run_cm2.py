import json
import logging
from os import path
from pathlib import Path
import time
from reconstruction import reconstruct
import fnmatch
from  pprint import pprint
from extension import reconstruct_slices, addition 

from futils import timeit
from tqdm import tqdm

from matching import Library, Sequence, match_library

# Logging configuration
current_file = path.basename(__file__).split(".")[0]

LYCOPENE = False

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
    extension_output_filename,
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

    # Get sequences TU in order of template
    with open(template_json_file) as json_file:
        template = json.load(json_file)

    # Get order of target sequences from the template slices order
    if LYCOPENE:
        slice_primer = {
                "revE": "CrtE",
                "revI": "CrtI",
                "revB": "CrtB",
            }
    else:
        slice_primer = {
            'TU-1': "tu1",
            'TU-2': "tu2",
            'TU-3': "tu3",
            'TU-4': "tu4",
            'TU-5': "tu5",
        }
    slices = template['template_slices']
    order_of_targets_paths = []
    order_of_targets_names = [] 
    for slice in slices:
        name = slice['name']
        primer = slice_primer[name]
        for f in seq_filenames:
            pattern = "*"+primer+"*"
            print(f, pattern)
            print("Name: ",f.name.split('.')[0])
            if f.match(pattern):
                order_of_targets_names.append(f.name.split('.')[0])
                order_of_targets_paths.append(f)

    print(order_of_targets_paths)
    print(order_of_targets_paths)
    # Loop over the sequences
    r = []
    for filename in order_of_targets_paths:
        sq = Sequence(filename)
        json_to_output = {}
        json_to_output["target"] = sq.name

        # Logging
        logging.info(f"Target sequence: {sq.name}")

        # Get libs from template
        libs = get_slices_libs(template)

        # Primer
        print(sq.name)
        #TODO name shit
        if LYCOPENE:
            primer = sq.name.split("_")[-2] # "BASIC_construct_UTR1-RBS-A12-UTR2-RBS-A12-UTR3-RBS-A12_CrtI_01"
            print("PRIMER:", primer)
            if primer == "CrtE":
                libs_to_match = libs['revE']
            elif primer == "CrtB":
                libs_to_match = libs['revB']
            elif primer == "CrtI":
                libs_to_match = libs['revI']
            else:
                raise OSError("Primer not found",sq.name)
        else:
            primer = sq.name.split("_")[2] # "vio_000_tu5"
            if primer == "tu1":
                libs_to_match = libs['TU-1']
            elif primer == "tu2":
                libs_to_match = libs['TU-2']
            elif primer == "tu3":
                libs_to_match = libs['TU-3']
            elif primer == "tu4":
                libs_to_match = libs['TU-4']
            elif primer == "tu5":
                libs_to_match = libs['TU-5']
            else:
                raise OSError("Primer not found",sq.name)

        # Match sequence
        matches = match_libs(sq, libs_to_match, threshold=threshold)
        json_to_output["matches"] = matches
        r.append(json_to_output)

    # Write output result in JSON file
    with open(match_output_filename, "w") as filename:
        json.dump(r, filename, indent=2, separators=(",", ":"))

    # Reconstruction
    reconstruction_result = reconstruct(r)

    # Write reconstruction result in JSON file
    with open(reconstruction_output_filename, "w") as filename:
        json.dump(reconstruction_result, filename, indent=2, separators=(",", ":"))

    # Extension and Addition to reconstruct the full pathway
    print(reconstruction_result)
    ex = reconstruct_slices(reconstruction_result, template, order_of_targets_names)
    print(ex)
    ex_res = {
            'target_slices': order_of_targets_names,
            'full_reconstruction': ex
            }
    with open(extension_output_filename , "w") as filename:
        json.dump(ex_res, filename, indent=2, separators=(",", ":"))


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

    msg = f"{test_name} - Targets:{targets['seq_dir_names']} - Template:{path.basename(candidates)} - Nb runs: {nbloop}"
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
        reconstruction_output_filename = path.join(
            OUTPUT_DIR, reconstruction_output_filename
        )

        # Final results filename
        extension_output_filename = f"{timestr}-extension-{current_file}-{test_id}-run-{i+1}-from-{nbloop}.json"
        extension_output_filename =  path.join(
            OUTPUT_DIR, extension_output_filename
        )


        # Iterate and match libs
        iter_all_seq(
            targets,
            candidates,
            match_output_filename,
            reconstruction_output_filename,
            extension_output_filename,
            threshold,
        )

# /////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////

def test_algo1_1_target_hard_th99():
    """
    Test Algo1
    1 Target Violacein (B0030, B0030, B0030, B0030, B0030) (hard)
    Candidate Template
    Threshold 0.99
    """

    test_params = {
        "name": "Algo1 - 1 Target (5*RBS B0030) - 1 against All - Threshold 0.99 ",
        "id": "vio-1-target-B0030-B0030-B0030-B0030-B0030-th99",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/violacein-basic-ass",
            "seq_dir_names": ["output/rbs_one_hard/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_violacein_02.json",
        "nbloop": 10,
        "threshold": 0.99,
    }

    run_test(test_params)


def test_algo1_1_target_easy_th75():
    """
    Test Algo1
    1 Target Violacein (B0030, B0031, B0032, B0033, B0064) (easy)
    Candidate Template
    Threshold 0.75
    """

    test_params = {
        "name": "Algo1 - 1 Target (RBS30,31,32,33,64) - 1 against all - Threshold 0.75 ",
        "id": "vio-1-target-B0030-B0031-B0032-B0033-B0064-th75",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/violacein-basic-ass",
            "seq_dir_names": ["output/rbs_one_easy/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_violacein_02.json",
        "nbloop": 10,
        "threshold": 0.75,
    }

    run_test(test_params)


def test_algo1_1_target_easy_th99():
    """
    Test Algo1
    1 Target RBS 30 31 32 33 64 (easy)
    Candidate Template
    Threshold 0.99
    """

    test_params = {
        "name": "Algo1 - 1 Target (RBS30,31,32,33,64) - 1 against all - Threshold 0.99 ",
        "id": "vio-1-target-B0030-B0031-B0032-B0033-B0064-th99",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/violacein-basic-ass",
            "seq_dir_names": ["output/rbs_one_easy/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_violacein_02.json",
        "nbloop": 10,
        "threshold": 0.99,
    }

    run_test(test_params)

# /////////////////////////////////////////////////////////////////////////////////////////

def test_algo1_1_target_1to1_hard_th75():
    """
    Test Algo1
    1 Target Violacein (B0030, B0030, B0030, B0030, B0030) (hard)
    Candidate Template
    Threshold 0.75
    """

    test_params = {
        "name": "Algo1 - 1 Target (5*RBS B0030) - 1 against 1 - Threshold 0.75",
        "id": "vio-1-target-B0030-B0030-B0030-B0030-B0030-1to1-th75",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/violacein-basic-ass",
            "seq_dir_names": ["output/rbs_one_hard/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_violacein_02_one_hard.json",
        "nbloop": 10,
        "threshold": 0.75,
    }

    run_test(test_params)


def test_algo1_1_target_1to1_hard_th99():
    """
    Test Algo1
    1 Target Violacein (B0030, B0030, B0030, B0030, B0030) (hard)
    Candidate Template
    Threshold 0.99
    """

    test_params = {
        "name": "Algo1 - 1 Target (5*RBS B0030) - 1 against 1 - Threshold 0.99",
        "id": "vio-1-target-B0030-B0030-B0030-B0030-B0030-1to1-th99",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/violacein-basic-ass",
            "seq_dir_names": ["output/rbs_one_hard/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_violacein_02_one_hard.json",
        "nbloop": 10,
        "threshold": 0.99,
    }

    run_test(test_params)


def test_algo1_1_target_1to1_easy_th75():
    """
    Test Algo1
    1 Target Violacein (B0030, B0031, B0032, B0033, B0064) (easy)
    Candidate Template
    Threshold 0.75
    """

    test_params = {
        "name": "Algo1 - 1 Target (RBS30,31,32,33,64) - 1 against 1 - Threshold 0.75",
        "id": "vio-1-target-B0030-B0031-B0032-B0033-B0064-1to1-th75",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/violacein-basic-ass",
            "seq_dir_names": ["output/rbs_one_easy/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_violacein_02_one_easy.json",
        "nbloop": 10,
        "threshold": 0.75,
    }

    run_test(test_params)


def test_algo1_1_target_1to1_easy_th99():
    """
    Test Algo1
    1 Target RBS 30 31 32 33 64 (easy)
    Candidate Template
    Threshold 0.99
    """

    test_params = {
        "name": "Algo1 - 1 Target (RBS30,31,32,33,64) - 1 against 1 - Threshold 0.99",
        "id": "vio-1-target-B0030-B0031-B0032-B0033-B0064-1to1-th99",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/violacein-basic-ass",
            "seq_dir_names": ["output/rbs_one_easy/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_violacein_02_one_easy.json",
        "nbloop": 10,
        "threshold": 0.99,
    }

    run_test(test_params)


# /////////////////////////////////////////////////////////////////////////////////////////

def test_algo1_2_targets_th99():
    """
    Test Algo1
    2 Targets
    Candidate Template
    Threshold 0.75
    Candidate Template
    """

    test_params = {
        "name": "Algo1 - 2 Targets - Candidate Template",
        "id": "2-targets-template-th99",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/violacein-basic-ass",
            "seq_dir_names": ["output/rbs_two/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_violacein_02.json",
        "nbloop": 10,
        "threshold": 0.99,
    }

    run_test(test_params)


def test_algo1_2_targets_th75():
    """
    Test Algo1
    2 Targets
    Candidate Template
    Threshold 0.75
    """

    test_params = {
        "name": "Algo1 - 2 Targets - Candidate Template",
        "id": "2-targets-template-th75",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/violacein-basic-ass",
            "seq_dir_names": ["output/rbs_two/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_violacein_02.json",
        "nbloop": 10,
        "threshold": 0.75,
    }

    run_test(test_params)


def test_algo1_10_targets():
    """
    Test Algo1
    10 Target
    Candidate Template
    """

    test_params = {
        "name": "Algo1 - 10 Targets - Candidate Template",
        "id": "target-10-template",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/violacein-basic-ass",
            "seq_dir_names": ["output/rbs_ten/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_violacein_02.json",
        "nbloop": 1,
    }

    run_test(test_params)


def main():
    """
    Main
    """
    pass


if __name__ == "__main__":
    main()
