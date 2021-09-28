from os import path
from futils import timeit, read_json
from pprint import pprint
import json

LYCOPENE = True

TEMPLATE_DIR = "/data/Imperial/src/matching/templates/"
TEMPLATE_FILE = "template_lycopene_sanger.json"
TEMPLATE = path.join(TEMPLATE_DIR, TEMPLATE_FILE)

RES_DIR = "/data/Imperial/src/matching/output_results/"

MAT_2 = "20210807-231259-matching-results-run_cm2_2-cm2-lycopene-1target-UTR1-RBS-A12-UTR2-RBS-A12-UTR3-RBS-A12-CrtI-th05-run-1-from-1.json"

REC_2 = "20210807-231259-reconstruction-run_cm2_2-cm2-lycopene-1target-UTR1-RBS-A12-UTR2-RBS-A12-UTR3-RBS-A12-CrtI-th05-run-1-from-1.json"
REC_2 = "20210812-111500-reconstruction-run_cm2_2-cm2-lycopene-100targets-sanger-th09-run-1-from-1.json"

MATCHES_2 = path.join(RES_DIR, MAT_2)
RECS_2 = path.join(RES_DIR, REC_2)


@timeit
def reconstruct_slices(rec, template, targets):
    """
    Extend

    targets = [
    {
        "name":"TargetABC",
        "slices": {
             "revE": "BASIC_construct_UTR1-RBS-A12-UTR2-RBS-A12-UTR3-RBS-A12_CrtE_01.seq",
             "revI": "BASIC_construct_UTR1-RBS-A12-UTR2-RBS-A12-UTR3-RBS-A12_CrtI_01.seq",
             "revB": "BASIC_construct_UTR1-RBS-A12-UTR2-RBS-A12-UTR3-RBS-A12_CrtB_01.seq"
        }
    },
    ...
    ]


    """
    template_size = len(template["template"]["structure"])
    template_slices = template["template_slices"]

    dd = {}
    l = []
    for slice in rec:
        print(slice)
        ll = []
        for ss in slice:
            d = {}
            sp = ss["path"]
            # Filter dict with wanted keys only 'name' and 'score'
            tr = [{k: d[k] for k in ["name", "score"]} for d in sp]
            ex = extend_slice(ss["slice"], template_slices, template_size, tr)
            slice_name = ss["slice_name"]
            d[slice_name] = ex
            d["score"] = ss["score"]
            d["name"] = slice_name
            ll.append(d)
        best = sorted(ll, key=lambda d: d["score"])
        # check if list is not empty
        if best:
            l.append(best[-1])
    print("Best:", l)
    r = addition(l, template_size)
    return r


def addition(es, template_size):
    """
    Addition
    """

    # TODO Change this, no hardcode
    target_json_file = (
        "/data/Imperial/src/matching/targets/target_lycopene_sanger_10.json"
    )

    with open(target_json_file) as json_file:
        target = json.load(json_file)

    targets = target["targets"]

    print("/" * 80)
    pprint(es)
    print("/" * 80)
    u = []
    for target in targets:
        li = []
        for k, v in target["slices"].items():
            print(k, v)
            v = path.basename(v)
            bn = v.split(".")[0]
            li.append(bn)
        print("Li:", li)
        addi = []
        for col in range(template_size):
            print("Part:", col + 1)
            l = []
            for k, v in enumerate(li):
                print("Slice:", k, v)
                a = []
                for i in es:
                    print("Name:", i["name"])
                    if v == i["name"]:
                        print("Match")
                        print(i[i["name"]])
                        a = i[i["name"]]
                if a:
                    l.append(a[col])
            print("L", l)
            cand = sorted(l, key=lambda k: k["score"])
            if cand:
                print(cand)
                print("-" * 80)
                addi.append(cand[-1])
        u.append({target["name"]: addi})
    return u


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

    # target_json_file = "/data/Imperial/src/matching/targets/target_lycopene_1.json"
    target_json_file = (
        "/data/Imperial/src/matching/targets/target_lycopene_sanger_A.json"
    )

    with open(target_json_file) as json_file:
        target = json.load(json_file)

    targets = target["targets"]

    r = reconstruct_slices(rec, template, targets)
    print("Construct: ", r)
    pprint(r)


if __name__ == "__main__":
    main()
