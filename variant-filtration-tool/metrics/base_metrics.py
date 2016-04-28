'''
Base class for CWL metrics tool
'''
import hashlib

from cdis_pipe_utils import time_util

class CWLMetricsTool(object):
    def __init__(self, time_file, normal_id, tumor_id, input_uuid, output_uuid, case_id, engine):
        self.time_file   = time_file
        self.normal_id   = normal_id
        self.tumor_id    = tumor_id
        self.input_uuid  = input_uuid
        self.output_uuid = output_uuid
        self.case_id     = case_id
        self.engine      = engine

    def get_time_metrics(self):
        time_str = None
        with open(self.time_file, 'rb') as fh:
            time_str = fh.read()
        return time_util.parse_time(time_str)

    def add_metrics(self):
        pass

class CWLMetricsMd5Tool(object):
    BLOCKSIZE = 65536

    def __init__(self, time_file, normal_id, tumor_id, input_uuid, output_uuid, case_id, engine, input_file):
        self.time_file   = time_file
        self.normal_id   = normal_id
        self.tumor_id    = tumor_id
        self.input_uuid  = input_uuid
        self.output_uuid = output_uuid
        self.case_id     = case_id
        self.engine      = engine
        self.input_file  = input_file

    def get_time_metrics(self):
        time_str = None
        with open(self.time_file, 'rb') as fh:
            time_str = fh.read()
        return time_util.parse_time(time_str)

    def add_metrics(self):
        pass

    def get_gz_md5(self):
        hasher = hashlib.md5()
        with open(self.input_file, 'rb') as afile:
            buf = afile.read(self.BLOCKSIZE)
            while len(buf) > 0:
                hasher.update(buf)
                buf = afile.read(self.BLOCKSIZE)
        return hasher.hexdigest()
