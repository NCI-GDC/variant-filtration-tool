"""Utility functions for the gdc_filtration_tools
package.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""

from typing import cast

from typing_extensions import Literal


def get_pysam_outmode(fname: str) -> Literal["r", "w", "wh", "rb", "wb", "wbu", "wb0"]:
    """
    Based on the filename returns wz etc.

    :param fname: the output filename
    :return: string pysam mode
    """
    mode = "w" if fname.endswith("gz") else "w"
    return cast(Literal["r", "w", "wh", "rb", "wb", "wbu", "wb0"], mode)
