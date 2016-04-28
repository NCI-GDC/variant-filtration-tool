'''
Metrics table class for fpfilter
'''
import gzip
from sqlalchemy import Column, Integer

from metrics.mixins import CustomToolMd5TypeMixin
from metrics.base_metrics import CWLMetricsMd5Tool

from cdis_pipe_utils import postgres


class FPFilterMetricsTable(CustomToolMd5TypeMixin, postgres.Base):
    total_variants = Column(Integer)

    __tablename__ = 'fpfilter_metrics'

class FPFilterMetricsTool(CWLMetricsMd5Tool):
    def __init__(self, time_file, normal_id, tumor_id, input_uuid, output_uuid, case_id, engine, input_file):
        super(FPFilterMetricsTool,self).__init__(time_file, normal_id, tumor_id, 
                                                     input_uuid, output_uuid, case_id, 
                                                     engine, input_file)
        self.tool  = 'fpfilter'
        self.files = [normal_id, tumor_id]

    def add_metrics(self):
        time_metrics = self.get_time_metrics()
        md5          = self.get_gz_md5()
        nvar         = get_variant_counts()
        metrics      = FPFilterMetricsTable(case_id      = self.case_id,
                                       vcf_id            = self.output_uuid,
                                       src_vcf_id        = self.input_uuid,
                                       tool              = self.tool,
                                       files             = self.files,
                                       systime           = time_metrics['system_time'],
                                       usertime          = time_metrics['user_time'],
                                       elapsed           = time_metrics['wall_clock'],
                                       cpu               = time_metrics['percent_of_cpu'],
                                       max_resident_time = time_metrics['maximum_resident_set_size'],
                                       md5               = md5,
                                       total_variants    = nvar)
        postgres.create_table(self.engine, metrics)
        postgres.add_metrics(self.engine, metrics)

    def get_variant_counts(self):
        reader = gzip.open(self.input_file, 'rt') if self.input_file.endswith('.gz') else open(self.input_file, 'r')
        
        total = 0
        for line in reader:
            if line.startswith('#'): continue
            else: total += 1
        return total 
