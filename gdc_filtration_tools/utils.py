"""Utility functions for the gdc_filtration_tools
package.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""


def get_pysam_outmode(fname: str):
    """
    Based on the filename returns wz etc.

    :param fname: the output filename
    :return: string pysam mode
    """
    mode = "wz" if fname.endswith("gz") else "w"
    return mode
