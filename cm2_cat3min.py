from run_cm2_cat import run_test


def test_cm2_cat3min():
    """
    Test CM_2 cat3 min
    Violacein cat 2x min with only RBS 30
    Threshold: 0.99
    """

    test_params = {
        "name": "CM2 cat3 min",
        "id": "cm2-cat3-min",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/cmatch/supplementary/Testing_Algorithm_CM_2/Viocat/",
            "seq_dir_names": ["targetcatx3/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_cm2_cat3min.json",
        "threshold": 0.99,
        "nbloop": 10,
    }

    run_test(test_params)


def main():
    """
    Main
    """
    test_cm2_cat3min()


if __name__ == "__main__":
    main()
