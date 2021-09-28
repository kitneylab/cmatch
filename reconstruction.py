import json
import itertools
from time import sleep, strftime
from os import path
from pprint import pprint
from futils import timeit, read_json

import logging
from statistics import geometric_mean


def compute_scores(paths):
    """
    Returns final score of the candidate pathway

    Arguments:
        paths list[str]

    Returns:
        list[float]  scores
    """
    scores = []
    for p in paths:
        gm = geometric_mean([e["score"] for e in p])
        scores.append(gm)
    return scores


def construct_names(paths):
    """ """
    names = []
    for p in paths:
        names.append("-".join([e["name"] for e in p]))
    return names


@timeit
def reconstruct(matches):
    """
    Reconstruction

    Args:
        matches (dict): JSON object

    Returns:
        list of dict containing the target, reconstruct, score and parts list
            ex: [ {
                    'target': 'vio-B0030-B0031-B0032-B0033-B0064',
                    'reconstruct': 'J23101-B0030-VioA-B0015-J23101-B0031-VioB-B0015-J23101-B0032-VioC-B0015-J23101-B0033-VioD-B0015-J23101-B0064-VioE-B0015',
                    'score': 20.0,
                    'path': [
                        {'name': 'J23101', 'score': 1.0, 'start': 4, 'length': 35, 'end': 39},
                        {'name': 'B0030', 'score': 1.0, 'start': 43, 'length': 15, 'end': 58},
                        {'name': 'VioA', 'score': 1.0, 'start': 62, 'length': 1293, 'end': 1355},
                        {'name': 'B0015', 'score': 1.0, 'start': 1359, 'length': 129, 'end': 1488},
                        {'name': 'J23101', 'score': 1.0, 'start': 1492, 'length': 35, 'end': 1527},
                        {'name': 'B0031', 'score': 1.0, 'start': 1531, 'length': 14, 'end': 1545},
                        {'name': 'VioB', 'score': 1.0, 'start': 1549, 'length': 3033, 'end': 4582},
                        {'name': 'B0015', 'score': 1.0, 'start': 4586, 'length': 129, 'end': 4715},
                        {'name': 'J23101', 'score': 1.0, 'start': 4719, 'length': 35, 'end': 4754},
                        {'name': 'B0032', 'score': 1.0, 'start': 4758, 'length': 13, 'end': 4771},
                        {'name': 'VioC', 'score': 1.0, 'start': 4775, 'length': 1326, 'end': 6101},
                        {'name': 'B0015', 'score': 1.0, 'start': 6105, 'length': 129, 'end': 6234},
                        {'name': 'J23101', 'score': 1.0, 'start': 6238, 'length': 35, 'end': 6273},
                        {'name': 'B0033', 'score': 1.0, 'start': 6277, 'length': 11, 'end': 6288},
                        {'name': 'VioD', 'score': 1.0, 'start': 6292, 'length': 1158, 'end': 7450},
                        {'name': 'B0015', 'score': 1.0, 'start': 7454, 'length': 129, 'end': 7583},
                        {'name': 'J23101', 'score': 1.0, 'start': 7587, 'length': 35, 'end': 7622},
                        {'name': 'B0064', 'score': 1.0, 'start': 7626, 'length': 12, 'end': 7638},
                        {'name': 'VioE', 'score': 1.0, 'start': 7642, 'length': 612, 'end': 8254},
                        {'name': 'B0015', 'score': 1.0, 'start': 8258, 'length': 129, 'end': 8387}
                        ]
                    }
                ]

    """
    # Read the JSON file with all the matches to reconstruct
    # targets = read_matches(matches)

    # Read the input list directly
    targets = matches
    target_reconstructions = []
    # Reconstruct each target
    for target in targets:
        print("Target:", target["target"])
        libs = target["matches"]
        candidates = []

        # Root
        paths = []
        for e in libs[0]["candidates"]:
            paths.append([e])

        print(strftime("%Y%m%d-%H%M%S"))
        print("Depth:", 0)
        print("\tnb paths:", len(paths))

        for i in range(1, len(libs), 1):
            # Add new lib
            np = []
            for pa in paths:
                aa = sorted(
                    libs[i]["candidates"], key=lambda d: d["score"], reverse=True
                )
                print(aa)
                aa = aa[0:1]
                # TODO verify highest score
                for e in aa:
                    new = pa.copy()
                    new.append(e)
                    np.append(new)
            print(np)
            # Prune
            paths = []

            print("Depth:", i)
            print("\tnb paths:", len(np))
            for p in np:
                print(p)
                if p[i - 1]["end"] <= p[i]["start"]:
                    print("ADDDDING:", p)
                    paths.append(p)
            print("\tafter pruning:", len(paths))
        scores = compute_scores(paths)
        names = construct_names(paths)
        r = []
        for i in range(len(paths)):
            d = {
                "target": target["target"],
                "reconstruct": names[i],
                "score": scores[i],
                "path": paths[i],
            }
            r.append(d)
        target_reconstructions.append(r)
    return target_reconstructions


def main():
    """
    Main
    """
    # Logging configuration
    current_file = path.basename(__file__).split(".")[0]
    logging.basicConfig(
        format="%(asctime)s:%(levelname)s: %(message)s",
        filename=f"logs/{current_file}.log",
        encoding="utf-8",
        level=logging.DEBUG,
    )

    # TEMPLATE = "/data/Imperial/src/matching/templates/template_violacein_02.json"
    # RES_DIR = "/data/Imperial/src/matching/output_results/"
    # RES = "matching-results-run_algo1-2-targets-template-run-0-20210804-130134.json"
    # MATCHES = path.join(RES_DIR, RES)
    TEMPLATE = "/data/Imperial/src/matching/templates/template_lycopene_sanger.json"
    RES_DIR = "/data/Imperial/src/matching/output_results/"
    RES = "20210807-185353-matching-results-run_cm2_2-cm2-lycopene-1target-UTR1-RBS-A12-UTR2-RBS-A12-UTR3-RBS-A12-CrtI-th05-run-1-from-1.json"
    MATCHES = path.join(RES_DIR, RES)

    r = reconstruct(read_json(MATCHES))
    print("Reconstruction:", r)
    print("Length:", len(r))


if __name__ == "__main__":
    main()
