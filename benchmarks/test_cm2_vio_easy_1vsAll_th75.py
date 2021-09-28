from run_cm2 import run_test


def test_cm2_vio_easy_1vsAll_th75():
    """
    Test CM_2
    Violacein
    1 vs All 
    Ease B0030-B0064
    Threshold: 0.75
    """

    test_params = {
        "name": "CM_2: Violacein - 1 vs All - Easy B0030-B0064 - Threshhold: 0.75",
        "id": "vio-1-easy-1vsAll-B0030-B0031-B0032-B0033-B0064-th75",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/violacein-basic-ass",
            "seq_dir_names": ["output/rbs_tu_easy/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_violacein_02_tu_1vsAll.json",
        "nbloop": 10,
        "threshold": 0.75,
    }

    run_test(test_params)


def main():
    """
    Main
    """
    test_cm2_vio_easy_1vsAll_th75()


if __name__ == "__main__":
    main()
