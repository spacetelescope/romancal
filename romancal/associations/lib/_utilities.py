"""General Utilities"""

import logging

from romancal.lib.basic_utils import parse_visitID as parse_visitID

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def mk_level3_asn_name(
    visit_id, output_file_root, filter_id, release_product, product_type, patch_name
):
    """Construct an association file name based on the visit_id and product"""

    parsed_visit_id = parse_visitID(visit_id)

    sep = "_"

    product_name_mapping = {
        "visit": "v"
        + parsed_visit_id["Execution"]
        + parsed_visit_id["Pass"]
        + parsed_visit_id["Segment"]
        + parsed_visit_id["Observation"]
        + parsed_visit_id["Visit"],
        "daily": "d"
        + parsed_visit_id["Execution"]
        + parsed_visit_id["Pass"]
        + parsed_visit_id["Segment"],
        "pass": "p" + parsed_visit_id["Execution"] + parsed_visit_id["Pass"],
        "full": "full",
        "user": "user",
    }

    pr_name = product_name_mapping.get(product_type, "unknown")

    asn_file_name = (
        output_file_root
        + sep
        + release_product
        + sep
        + pr_name
        + sep
        + patch_name
        + sep
        + filter_id
    )

    return asn_file_name
