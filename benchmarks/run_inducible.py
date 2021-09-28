from os import path, listdir
from matching import *
from Bio import SeqIO, Seq
from pprint import pprint
from pathlib import Path

from templates import template_inducible as template_inducible


def match_libs(seq, libs, threshold=0.5):
    """
    Match libs
    with the sequence from probab trace
    """
    for lib in libs:
        print("-" * 80)
        print(f"{lib.name}:")
        candidates = match_library(seq, lib, threshold)
        for l in candidates:
            for c in l:
                s = f"{c.name} {c.score} {c.start} {c.length} {c.end}"
                print(s)


def match_libs_proba(seq, libs, threshold=0.5):
    """
    Match libs
    with the sequence from probab trace
    """
    for lib in libs:
        print("-" * 80)
        print(f"{lib.name}:")
        candidates = match_library_proba(seq, lib, threshold)
        for l in candidates:
            for c in l:
                s = f"{c.name} {c.score} {c.start} {c.length} {c.end}"
                print(s)


def get_sequences():
    """
    Return list of sequences
    """

    SEQUENCES_EXTENSION = ".seq"
    SEQUENCES_PATH = "/data/Imperial/repo-inducible/sequencing/"
    #SEQUENCES_PATH = "/data/Imperial/src/lyc-basic-ass-ind/"
    sequences = []
    seq_dir_names = [
        "seq-sourcebio-21-462006601-001",
        "seq-sourcebio-25-462779601-001",
        "seq-sourcebio-28-462779601-47f67342-001",
        "seq-sourcebio-39-462743501-001",
    ]
#    seq_dir_names = [ "output" ]
    for seq_dir in seq_dir_names:
        seqs = Path(path.join(SEQUENCES_PATH, seq_dir)).rglob(
            "*{0}".format(SEQUENCES_EXTENSION)
        )
        sequences.append(seqs)
    return sequences


def get_sub_templates_libs():
    """
    Get sub-templates libraries
    """

    libs_E = []
    for pos in template_inducible.sub_E:
        lib = template_inducible.template["structure"][pos - 1]["library_source"]
        libs_E.append(Library(eval("template_inducible." + lib)))
    libs_I = []
    for pos in template_inducible.sub_I:
        lib = template_inducible.template["structure"][pos - 1]["library_source"]
        libs_I.append(Library(eval("template_inducible." + lib)))
    libs_B = []
    for pos in template_inducible.sub_B:
        lib = template_inducible.template["structure"][pos - 1]["library_source"]
        libs_B.append(Library(eval("template_inducible." + lib)))
    return libs_E, libs_I, libs_B


def iter_all_seq_proba():
    """
    Iterate over inducible sequences
    """

    # Get sequences to match
    sequences = get_sequences()

    # Get sub-templates libraries for primers E, I and B
    libs_E, libs_I, libs_B = get_sub_templates_libs()

    # Get the filenames in a list and not this freakin generator
    seq_filenames = []
    for seq in sequences:
        for filename in seq:
            seq_filenames.append(filename)

    # Loop over the sequences
    for filename in sorted(seq_filenames):
        sq = Sequence(filename)
        print("-" * 80)
        print("-" * 80)
        print(sq.name)
        # Get libraries to match according to the primer
        #primer = sq.name.split("_")[2]
        primer = sq.name.split("_")[3]
        print(primer)
        if primer == "CrtE":
            libs_to_match = libs_E
        elif primer == "CrtI":
            libs_to_match = libs_I
        elif primer == "CrtB":
            libs_to_match = libs_B
        else:
            raise OSError("Primer not found",sq.name)
        # Match sequence
        # match_libs_proba(sq, libs_to_match, threshold=0.1)
        match_libs(sq, libs_to_match, threshold=0.1)


def test_i01():
    """
    Test for construct i01
    """

    SEQUENCES_PATH = (
        "/data/Imperial/repo-inducible/sequencing/seq-sourcebio-21-462006601-001"
    )
    sequences = []
    seq_names = [
        "462006601_i01_CrtE_B04.seq",
        "462006601_i01_CrtI_C04.seq",
        "462006601_i01_CrtB_A04.seq",
    ]

    for seq_name in seq_names:
        f = path.join(SEQUENCES_PATH, seq_name)
        sequences.append(f)

    # Libs to match
    lib_prom = Library(
        template_inducible.library_inducible_parts_promoters
    )  # promoters
    lib_cds_E = Library(template_inducible.library_inducible_parts_cds_E)  # crtE
    lib_cds_I = Library(template_inducible.library_inducible_parts_cds_I)  # crtI
    lib_cds_B = Library(template_inducible.library_inducible_parts_cds_B)  # crtB
    lib_utr1_rbs = Library(template_inducible.library_inducible_parts_utr1_rbs)
    lib_utr2_rbs = Library(template_inducible.library_inducible_parts_utr2_rbs)
    lib_utr3_rbs = Library(template_inducible.library_inducible_parts_utr3_rbs)

    # Get sub-templates libraries for primers E, I and B
    libs_E = []
    for pos in template_inducible.sub_E:
        lib = template_inducible.template["structure"][pos - 1]["library_source"]
        libs_E.append(Library(eval("template_inducible." + lib)))
    libs_I = []
    for pos in template_inducible.sub_I:
        lib = template_inducible.template["structure"][pos - 1]["library_source"]
        libs_I.append(Library(eval("template_inducible." + lib)))
    libs_B = []
    for pos in template_inducible.sub_B:
        lib = template_inducible.template["structure"][pos - 1]["library_source"]
        libs_B.append(Library(eval("template_inducible." + lib)))

    # Loop over the sequences
    for filename in sequences:
        sq = Sequence(filename)
        print("-" * 80)
        print("-" * 80)
        print(sq.name)
        # Get libraries to match according to the primer
        primer = sq.name.split("_")[2]
        if primer == "CrtE":
            libs_to_match = libs_E
        elif primer == "CrtI":
            libs_to_match = libs_I
        elif primer == "CrtB":
            libs_to_match = libs_B
        # Match sequence
        # match_libs_proba(sq, libs_to_match, threshold=0.1)
        match_libs(sq, libs_to_match, threshold=0.1)


def get_all_sequences():
    SEQUENCES_EXTENSION = ".seq"
    SEQUENCES_PATH = "/data/Imperial/repo-inducible/sequencing/"
    sequences = []
    seq_dir_names = listdir(SEQUENCES_PATH)
    seq_dir_names.remove("downloads")  # Remove the downloads/ directory
    for seq_dir in seq_dir_names:
        seqs = Path(path.join(SEQUENCES_PATH, seq_dir)).rglob(f"*{SEQUENCES_EXTENSION}")
        sequences.append(seqs)
    num = []
    for s in sequences:
        for ss in s:
            bn = path.basename(str(ss))
            num.append(path.splitext(bn)[0])
    fn = list(set(num))
    print(fn, len(fn))
    sn = [f.split("_")[1] for f in fn]
    x = sorted(list(set(sn)))
    print(x, len(x))


def main():
    """
    Main
    """
    # test_i01()
    # get_all_sequences()
    iter_all_seq_proba()


if __name__ == "__main__":
    main()
