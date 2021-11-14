import json
import logging
from os import path
from pathlib import Path
import time
from reconstruction2 import reconstruct
import fnmatch
from pprint import pprint
from extension100 import reconstruct_slices, addition

from futils import timeit
from tqdm import tqdm

from matching import Library, Sequence, match_library

# Logging configuration
current_file = path.basename(__file__).split(".")[0]

LYCOPENE = True


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
    targets_param,
    template_json_file,
    match_output_filename,
    reconstruction_output_filename,
    extension_output_filename,
    threshold=0.99,
):
    """
    Iterate over target sequences
    """

    target_json_file = targets_param["file"]

    # Get targets from JSON file
    with open(target_json_file) as json_file:
        target = json.load(json_file)

    # Get sequences TU in order of template
    with open(template_json_file) as json_file:
        template = json.load(json_file)

    print("Target:", target)
    print("Template:", template)
    targets = target["targets"]
    slices = template["template_slices"]
    print("Slices:", slices)
    print("Targets_param:", targets_param)

    # Loop over the targets
    r = []
    order_of_targets_names = []
    for target in targets:
        name = target["name"]
        order_of_targets_names.append(name)
        # Logging
        logging.info(f"Target sequence: {name}")
        print("Target name:", name)

        # Match the slices
        for slice in slices:
            name = slice["name"]
            print("Slice name:", name)
            slice_file = target["slices"][name]
            print("Slice file:", slice_file)
            sq = Sequence(slice_file)
            json_to_output = {}
            print("Sq name:", sq.name)
            json_to_output["target"] = sq.name
            json_to_output["slice"] = name  # "revE"
            # Get libs from template
            libs = get_slices_libs(template)
            print("Libs:", libs)

            # Match slice
            libs_to_match = libs[name]
            matches = match_libs(sq, libs_to_match, threshold=threshold)
            json_to_output["matches"] = matches
            r.append(json_to_output)

            # Write output result in JSON file
            with open(match_output_filename, "w") as filename:
                json.dump(r, filename, indent=2, separators=(",", ":"))

    # Reconstruction
    reconstruction_result = reconstruct(r)
    print("Reconstruction:", reconstruction_result)

    # Write reconstruction result in JSON file
    with open(reconstruction_output_filename, "w") as filename:
        json.dump(reconstruction_result, filename, indent=2, separators=(",", ":"))

    # Extension and Addition to reconstruct the full pathway
    ex = reconstruct_slices(reconstruction_result, template, targets)
    print("Extension:", ex)
    with open(extension_output_filename, "w") as filename:
        json.dump(ex, filename, indent=2, separators=(",", ":"))


def run_test(test_params):
    """
    Run tests
    """

    targets = test_params["targets"]
    template = test_params["template"]
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

    msg = f"{test_name} - Targets:{path.basename(targets['file'])} - Template:{path.basename(template)} - Nb runs: {nbloop}"
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
        extension_output_filename = (
            f"{timestr}-extension-{current_file}-{test_id}-run-{i+1}-from-{nbloop}.json"
        )
        extension_output_filename = path.join(OUTPUT_DIR, extension_output_filename)

        # Iterate and match libs
        iter_all_seq(
            targets,
            template,
            match_output_filename,
            reconstruction_output_filename,
            extension_output_filename,
            threshold,
        )


def main():
    """
    Main
    """
    pass


if __name__ == "__main__":
    main()
