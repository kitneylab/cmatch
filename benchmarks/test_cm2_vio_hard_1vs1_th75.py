from run_cm2 import run_test


def test_cm2_vio_hard_1vs1_th75():
    """
    Test CM_2
    Violacein
    1 vs 1
    1 Target: Hard (B0030, B0030, B0030, B0030, B0030) (hard)
    Threshold: 0.75
    """

    test_params = {
        "name": "CM_2: Violacein - 1 vs 1 - Hard B0030..B0030 - Threshhold: 0.75",
        "id": "vio-1-target-B0030-B0030-B0030-B0030-B0030-th75",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/violacein-basic-ass",
            "seq_dir_names": ["output/rbs_tu_hard/"],
        },

        "candidates": "/data/Imperial/src/matching/templates/template_violacein_02_tu_1vs1_hard.json",
        "nbloop": 10,
        "threshold": 0.75,
    }

    run_test(test_params)


def main():
    """
    Main
    """
    test_cm2_vio_hard_1vs1_th75()


if __name__ == "__main__":
    main()
