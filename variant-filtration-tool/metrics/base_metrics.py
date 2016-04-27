'''
Base class for CWL metrics tool
'''
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
