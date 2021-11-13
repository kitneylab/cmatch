from run_cm2_cat import run_test


def test_cm2_cat2():
    """
    Test CM_2 cat2
    Violacein cat 2x
    Threshold: 0.99
    """

    test_params = {
        "name": "CM2 Vio 0000 cat2",
        "id": "cm2-vio-0000-cat2",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/cmatch/supplementary/Testing_Algorithm_CM_2/Viocat/",
            "seq_dir_names": ["targetcatx2/"],
        },
        "candidates": "/data/Imperial/src/matching/templates/template_cm2_cat2.json",
        "threshold": 0.99,
        "nbloop": 1,
    }

    run_test(test_params)


def main():
    """
    Main
    """
    test_cm2_cat2()


if __name__ == "__main__":
    main()
