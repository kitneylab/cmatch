from os import path
from futils import timeit, read_json
from pprint import pprint

LYCOPENE = True

TEMPLATE_DIR = "/data/Imperial/src/matching/templates/"
TEMPLATE_FILE = "template_violacein_02_tu_1vs1_easy.json"
TEMPLATE = path.join(TEMPLATE_DIR, TEMPLATE_FILE)

RES_DIR = "/data/Imperial/src/matching/output_results/"

# EASY
MAT_1 = (
    "20210805-151049-matching-results-run_algo2_rbs-1-target-template-run-1-from-1.json"
)
REC_1 = (
    "20210805-151049-reconstruction-run_algo2_rbs-1-target-template-run-1-from-1.json"
)
MATCHES_1 = path.join(RES_DIR, MAT_1)
RECS_1 = path.join(RES_DIR, REC_1)

# HARD
MAT_2 = "20210805-204000-matching-results-run_cm2-vio-1-target-B0030-B0030-B0030-B0030-B0030-th75-run-3-from-10.json"
REC_2 = "20210805-204000-reconstruction-run_cm2-vio-1-target-B0030-B0030-B0030-B0030-B0030-th75-run-3-from-10.json"
MATCHES_2 = path.join(RES_DIR, MAT_2)
RECS_2 = path.join(RES_DIR, REC_2)

TARGET_ORDER_1 = [
    "vio_A-0030_tu1",
    "vio_B-0031_tu2",
    "vio_C-0032_tu3",
    "vio_D-0033_tu4",
    "vio_E-0064_tu5",
]

TARGET_ORDER_2 = [
    "vio_A-0030_tu1",
    "vio_B-0030_tu2",
    "vio_C-0030_tu3",
    "vio_D-0030_tu4",
    "vio_E-0030_tu5",
]


def get_targets_order():
    """
    Returns target order according to the template
    """
    order = [
        "vio_A-0030_tu1",
        "vio_B-0031_tu2",
        "vio_C-0032_tu3",
        "vio_D-0033_tu4",
        "vio_E-0064_tu5",
    ]

    return order


def primer_from_slice(slice):
    """
    Return primer name from slice name
    """
    slice_primer = {
        "TU-1": "tu1",
        "TU-2": "tu2",
        "TU-3": "tu3",
        "TU-4": "tu4",
        "TU-5": "tu5",
    }
    return slice_primer[slice]


def slice_from_primer(primer):
    """
    Return slice name from primer name
    """

    if LYCOPENE:
        primer_slice = {
            "CrtE": "revE",
            "CrtI": "revI",
            "CrtB": "revB",
        }
    else:
        primer_slice = {
            "tu1": "TU-1",
            "tu2": "TU-2",
            "tu3": "TU-3",
            "tu4": "TU-4",
            "tu5": "TU-5",
        }
    return primer_slice[primer]


@timeit
def reconstruct_slices(rec, template, targets_order):
    """
    Extend
    """
    template_size = len(template["template"]["structure"])
    template_slices = template["template_slices"]
    # TODO change this shit below
    # order_of_targets = get_targets_order()
    order_of_targets = targets_order

    d = {}
    for ti, target_name in enumerate(order_of_targets):  # [ "vio_D-0033_tu4", .. ]
        # TODO change line below with a proper way
        if LYCOPENE:
            print("XXXXX", target_name)
            tu = target_name.split("_")[-2]
            print("YYYYY", tu)

        else:
            tu = target_name[-3:]
        slice = slice_from_primer(tu)
        for target in rec:
            # TODO
            print("Target: ", target, target_name)
            print("Target: ", target[0]["target"], target_name)
            if target[0]["target"] == target_name:
                target_result = target[0]["path"]
                # Filter dict with wanted keys only
                # dict_you_want = { your_key: old_dict[your_key] for your_key in your_keys }
                tr = [{k: d[k] for k in ["name", "score"]} for d in target_result]
                ex = extend_slice(slice, template_slices, template_size, tr)
                d[target_name] = ex
    print("Extended slices: ", d)
    r = addition(d, template_size)
    return r


def addition(extended_slices, template_size):
    """
    Addition
    """
    # print(extended_slices)
    addi = []
    for col in range(template_size):
        print("Part:", col + 1)
        l = []
        for k, v in extended_slices.items():
            print(k, v[col])
            l.append(v[col])
        cand = sorted(l, key=lambda k: k["score"])[-1]
        print(cand)
        print("-" * 80)
        addi.append(cand)
    return addi


def extend_slice(slice_name, slices, template_size, target_result):
    """
    Extend slice
    """
    ext_obj = {"name": "Empty", "score": 0.0}
    extended_list = [ext_obj] * template_size
    el = extended_list
    for slice in slices:
        if slice["name"] == slice_name:
            print("Slice:", slice_name)
            for i, n in enumerate(
                slice["template_slice"]
            ):  #  "template_slice":[1, 2, 3, 4]
                pos = n - 1
                el = el[:pos] + [target_result[i]] + el[pos + 1 :]
    return el


def main():
    """
    Main
    """

    rec = read_json(RECS_2)
    template = read_json(TEMPLATE)
    r = reconstruct_slices(rec, template, TARGET_ORDER_2)
    print("Construct: ", r)


if __name__ == "__main__":
    main()
