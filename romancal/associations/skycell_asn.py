"""Create an association based on skycells"""
import argparse
import sys

import json
import logging
import numpy as np

from romancal.associations import asn_from_list
import romancal.patch_match.patch_match as pm
import roman_datamodels as rdm


__all__ = ['skycell_asn']

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.setLevel('INFO')

def skycell_asn(self):
    """Create the associaton from the list"""
    all_patches = []
    file_list = []
    for file_name in self.parsed.filelist:
        cal_file = rdm.open(file_name)
        filter_id = cal_file.meta.instrument.optical_element
        file_patch_list = pm.find_patch_matches(cal_file.meta.wcs)
        logger.info(f"Patch List:   {file_name}  {file_patch_list[0]}")
        file_list.append([file_name, file_patch_list[0]])
        all_patches.append(file_patch_list[0])

    unique_patches = np.unique(np.concatenate(all_patches))
    for item in unique_patches:
        member_list = []
        patch_name = pm.PATCH_TABLE[item]['name']
        for a in file_list:
            if np.isin(item, a[1]):
                member_list.append(a[0])

        # grab all the wcs parameters needed for generate_tan_wcs
        projcell_info = dict([('pixel_scale', float(pm.PATCH_TABLE[item]['pixel_scale'])),
                              ('ra_cent' ,float(pm.PATCH_TABLE[item]['ra_projection_center'])),
                              ('dec_cent' ,float(pm.PATCH_TABLE[item]['dec_projection_center'])),
                              ('shiftx' , float(pm.PATCH_TABLE[item]['x0_projection'])),
                              ('shifty' , float(pm.PATCH_TABLE[item]['y0_projection'])),
                              ('nx', int(pm.PATCH_TABLE[item]['nx'])),
                              ('ny', int(pm.PATCH_TABLE[item]['ny'])),
                              ])
        program_id = member_list[0][1:6]
        root_asn_name = self.parsed.output_file_root
        product_type = self.parsed.product_type
        product_release = self.parsed.release_product
        suffix = "i2d"
        sep = "_"
        asn_file_name = root_asn_name + sep + patch_name + sep + product_type + \
            sep + filter_id + sep + product_release + sep + suffix
        with open(asn_file_name+"_asn.json", "w") as outfile:
            prompt_product_asn = asn_from_list.asn_from_list(member_list,
                                                             product_name = asn_file_name )
            prompt_product_asn["asn_type"] = "image"
            prompt_product_asn["program"] = program_id
            prompt_product_asn["target"] = patch_name
            #prompt_product_asn["skycell_wcs_info"] = projcell_info
            prompt_product_asn["skycell_wcs_info"] = json.dumps(projcell_info)
            _, serialized = prompt_product_asn.dump(format="json")
            outfile.write(serialized)


class Main():
    """Command-line interface for list_to_asn

    Parameters
    ----------
    args: [str, ...], or None
        The command line arguments. Can be one of
            - `None`: `sys.argv` is then used.
            - `[str, ...]`: A list of strings which create the command line
              with the similar structure as `sys.argv`
    """
    def __init__(self, args=None):
        if args is None:
            args = sys.argv[1:]
        if isinstance(args, str):
            args = args.split(" ")

        parser = argparse.ArgumentParser(
            description='Create an association from a list of files',
            usage='skycell_asn --product-type visit --release-product prompt *_cal.asdf -o r512',
        )

        parser.add_argument(
            '-o', '--output-file-root',
            type=str,
            required=True,
            help='Root string for file to write association to'
        )

        parser.add_argument(
            '-f', '--format',
            type=str,
            default='json',
            help='Format of the association files. Default: "%(default)s"'
        )

        parser.add_argument(
            '--product-type',
            type=str,
            default='visit',
            help='The product type when creating the association'
        )

        parser.add_argument(
            '--release-product',
            type=str,
            default='prompt',
            help='The release product when creating the association'
        )

        parser.add_argument(
            '-r', '--rule',
            type=str,
            default='DMS_ELPP_Base',
            help=(
                'The rule to base the association structure on.'
                ' Default: "%(default)s"'
            )
        )
        parser.add_argument(
            '-i', '--id',
            type=str,
            default='o999',
            help='The association candidate id to use. Default: "%(default)s"',
            dest='acid'
        )

        parser.add_argument(
            'filelist',
            type=str,
            nargs='+',
            help='File list to include in the association'
        )

        self.parsed = parser.parse_args(args=args)
        logger.info("Command-line arguments: %s", self.parsed)

        skycell_asn(self)

if __name__ == '__main__':

    Main()
