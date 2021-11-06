from run_cm2_10 import run_test

TEST_DIR="./supplementary/Real_Life_Example_Lycopene_Operon/"


def test_cm2_lycopene_10():
    """
    Test CM_2
    Lycopene Sanger
    10 targets
    Threshold: 0.7
    """

    test_params = {
        "name": "CM_2: Lycopene - 10 targets - Sanger - Threshold: 0.y",
        "id": "cm2-lycopene-10targets-sanger-th07",
        "targets": {
            "file": f"{TEST_DIR}/targets/target_lycopene_sanger_10.json",
        },
        "template": f"{TEST_DIR}/templates/template_lycopene_sanger.json",
        "nbloop": 1,
        "threshold": 0.7,
    }

    run_test(test_params)


def main():
    """
    Main
    """
    test_cm2_lycopene_10()


if __name__ == "__main__":
    main()
