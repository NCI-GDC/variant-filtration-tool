"""Utility functions for the gdc_filtration_tools
package.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""

from typing import Literal


def get_pysam_outmode(fname: str) -> Literal["w"]:
    """
    Based on the filename returns wz etc.

    :param fname: the output filename
    :return: string pysam mode
    """
    mode = "w" if fname.endswith("gz") else "w"
    return mode
