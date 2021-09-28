from run_cm2_2 import run_test


def test_cm2_lycopene_1():
    """
    Test CM_2
    Lycopene
    One target
    Threshold: 0.5
    """

    test_params = {
        "name": "CM_2: Lycopene - 1 target - UTR1-RBS-A12-UTR2-RBS-A12-UTR3-RBS-A12-CrtI - Threshold: 0.5",
        "id": "cm2-lycopene-1target-UTR1-RBS-A12-UTR2-RBS-A12-UTR3-RBS-A12-CrtI-th05",
        "targets": {
            "extension": ".seq",
            "sequences_path": "/data/Imperial/src/lyc-basic-ass-ind/",
            "seq_dir_names": ["output/primer/"],
            "file": "/data/Imperial/src/matching/targets/target_lycopene_1.json",
        },
        "candidates": "/data/Imperial/src/matching/templates/template_lycopene_sanger.json",
        "nbloop": 1,
        "threshold": 0.5,
    }

    run_test(test_params)


def main():
    """
    Main
    """
    test_cm2_lycopene_1()


if __name__ == "__main__":
    main()
