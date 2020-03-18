"""Testing utility functions."""
import os
import sys

from io import StringIO
from contextlib import contextmanager

from gdc_filtration_tools.logger import Logger


@contextmanager
def captured_output():
    """Captures stderr and stdout and returns them"""
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        Logger.setup_root_logger()
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


def get_test_data_path(fname: str):
    """
    Helper function to get full path to a test
    dataset give the basename.
    :param fname: test dataset file basename
    """
    path = os.path.join(os.path.dirname(__file__), "data/{0}".format(fname))
    return path


def cleanup_files(files):
    """
    Takes a file or a list of files and removes them.
    """
    def _do_remove(fil):
        if os.path.exists(fil):
            os.remove(fil)

    flist = []
    if isinstance(files, list):
        flist = files[:]
    else:
        flist = [files]

    for fil in flist:
        _do_remove(fil)
