from run_cm2 import run_test

TEST_DIR="./supplementary/Testing_Algorithm_CM_2/First_Test/"

def test_cm2_vio_easy_1vs1_th99():
    """
    Test CM_2
    Violacein
    1 vs 1
    Ease B0030-B0064
    Threshold: 0.99
    """

    test_params = {
        "name": "CM_2: Violacein - 1 vs 1 - Easy B0030-B0064 - Threshhold: 0.99",
        "id": "vio-1-easy-B0030-B0031-B0032-B0033-B0064-th99",
        "targets": {
            "extension": ".seq",
            "sequences_path": f"{TEST_DIR}/sequences", 
            "seq_dir_names": ["rbs_tu_easy/"],
        },
        "candidates": f"{TEST_DIR}templates/template_violacein_02_tu_1vs1_easy.json",
        "nbloop": 10,
        "threshold": 0.99,
    }

    run_test(test_params)


def main():
    """
    Main
    """
    test_cm2_vio_easy_1vs1_th99()


if __name__ == "__main__":
    main()
