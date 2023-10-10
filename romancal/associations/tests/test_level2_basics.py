"""Test basic usage of Level2 associations"""

from romancal.associations import generate, load_asn
from romancal.associations.main import Main
from romancal.associations.tests.helpers import (
    combine_pools,
    registry_level2_only,
    t_path,
)

# REGEX_LEVEL2 = r'(?P<path>.+)(?P<type>_rate(ints)?)(?P<extension>\..+)'
# REGEX_LEVEL2 = r'(?P<path>.+)(.*\_.*\_.*\_.*\_.*)'
REGEX_LEVEL2 = r"(?P<path>.+)(?P<type>_cal?)(?P<extension>\..+)"


def from_level2_schema():
    with open(t_path("data/asn_level2.json")) as asn_file:
        asn = load_asn(asn_file)
    return [asn]


def generate_from_pool(pool_path):
    """Generate associations from pools"""
    rules = registry_level2_only()
    pool = combine_pools(t_path(pool_path))
    asns = generate(pool, rules)
    return asns


def cmd_from_pool(pool_path, args):
    """Run commandline on pool

    Parameters
    ---------
    pool_path: str
        The pool to run on.

    args: [arg(, ...)]
        Additional command line arguments in the form `sys.argv`
    """
    full_args = [
        "--dry-run",
        "-D",
        "-r",
        t_path("../lib/rules_level2.py"),
        "--ignore-default",
    ]
    full_args.extend(args)
    result = Main(full_args, pool=pool_path)
    return result


def test_level2_productname():
    asns = generate_from_pool("data/pool_002_wfi_image.csv")
    for asn in asns:
        for product in asn["products"]:
            science = [
                member
                for member in product["members"]
                if member["exptype"] == "science" or member["exptype"] == "wfi_image"
            ]
            if asn["asn_rule"] == "Asn_Lv2Image":
                assert len(science) == 2
            if asn["asn_rule"] == "Asn_Lv2FOV":
                assert len(science) == 18


#            match = re.match(REGEX_LEVEL2, science[0]['expname'])
#            pdb.set_trace()
#            assert match.groupdict()['path'] == product['name']
