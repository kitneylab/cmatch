from run_cm2_cat import run_test


def test_cm2_cat3():
    """
    Test CM_2 cat3
    Violacein cat 3x
    Threshold: 0.99
    """

    test_params = {
        "name": "CM2 cat3",
        "id": "cm2-cat3",
        "targets": {
            "extension": ".seq",
            "sequences_path": "supplementary/Testing_Algorithm_CM_2/Viocat/",
            "seq_dir_names": ["targetcatx3/"],
        },
        "candidates": "supplementary/Testing_Algorithm_CM_2/Viocat/templates/template_cm2_cat3.json",
        "threshold": 0.99,
        "nbloop": 10,
    }

    run_test(test_params)


def main():
    """
    Main
    """
    test_cm2_cat3()


if __name__ == "__main__":
    main()
