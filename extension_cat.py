from os import path
from futils import timeit, read_json
from pprint import pprint

LYCOPENE = False


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
            "tu01": "TU-01",
            "tu02": "TU-02",
            "tu03": "TU-03",
            "tu04": "TU-04",
            "tu05": "TU-05",
            "tu06": "TU-06",
            "tu07": "TU-07",
            "tu08": "TU-08",
            "tu09": "TU-09",
            "tu10": "TU-10",
            "tu11": "TU-11",
            "tu12": "TU-12",
            "tu13": "TU-13",
            "tu14": "TU-14",
            "tu15": "TU-15",
        }
    return primer_slice[primer]


@timeit
def reconstruct_slices(rec, template, order_of_targets):
    """
    Extend
    """
    template_size = len(template["template"]["structure"])
    template_slices = template["template_slices"]
    print("Templat slices", template_slices)

    d = {}
    for ti, target_name in enumerate(order_of_targets):  # [ "vio_D-0033_tu4", .. ]
        # TODO change line below with a proper way
        if LYCOPENE:
            print("XXXXX", target_name)
            tu = target_name.split("_")[-2]
            print("YYYYY", tu)
        else:
            tu = target_name[-4:]  ## assumption on name (-4 because tu01.seq)
            print("Target name:", target_name, tu)
        slice = slice_from_primer(tu)
        for target in rec:
            # TODO
            # print("Target: ",target, target_name)
            # print("Target: ",target[0]["target"], target_name)
            if target[0]["target"] == target_name:
                target_result = target[0]["path"]
                print("MATCH")
                print("Target: ", target, target_name)
                print("Target: ", target[0]["target"], target_name)
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
    print("SliceName:", slice_name)
    print("Slices:", slices)
    ext_obj = {"name": "Empty", "score": 0.0}
    extended_list = [ext_obj] * template_size
    el = extended_list
    for slice in slices:
        print(slice["name"], slice_name)
        if slice["name"] == slice_name:
            print("X" * 80)
            print("X" * 80)
            print("Slice:", slice_name)
            for i, n in enumerate(
                slice["template_slice"]
            ):  #  "template_slice":[1, 2, 3, 4]
                pos = (n % template_size) - 1
                el = el[:pos] + [target_result[i]] + el[pos + 1 :]
    print("Extended Slice:", el)
    return el
