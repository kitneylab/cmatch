from run_cm2 import run_test


def test_cm2_ten():
    """
    Test CM_2  
    Violacein 
    Threshold: 0.99
    """
    TEST_DIR = "./supplementary/Testing_Algorithm_CM_2/"
    test_params = {
        "name": "CM2 Vio tenx",
        "id": "cm2-vio-tenx",
        "targets": {
            "extension": ".seq",
#            "sequences_path": f"{TEST_DIR}",
            "sequences_path": f"/data/Imperial/src/violacein-basic-ass/output/tudir",
            "seq_dir_names": [
                "0031-0032-0032-0033-0030",
                "0031-0032-0032-0064-0031",
                "0031-0032-0033-0032-0031",
                "0033-0030-0064-0064-0031",
                "0064-0064-0031-0064-0031",
                "0064-0064-0032-0031-0031",
                "0033-0031-0032-0033-0064",
                "0031-0033-0031-0032-0030",
                "0064-0064-0064-0064-0030",
                "0032-0064-0031-0064-0032",
               ],
        },
        "candidates": f"{TEST_DIR}/First_Test/templates/template_vio_tu_all.json",
        "threshold": 0.99,
        "nbloop": 1,
    }

    run_test(test_params)


def main():
    """
    Main
    """
    test_cm2_ten()


if __name__ == "__main__":
    main()
