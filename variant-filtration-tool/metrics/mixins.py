'''
Custom POSTGRES mixin classes
'''
from sqlalchemy import Column, Integer, String, Float
from sqlalchemy.dialects.postgresql import ARRAY

class CustomToolTypeMixin(object):
    ''' Gather timing metrics with input/output uuids '''
    id                = Column(Integer, primary_key=True)
    case_id           = Column(String)
    vcf_id            = Column(String)
    src_vcf_id        = Column(String)
    tool              = Column(String)
    files             = Column(ARRAY(String))
    systime           = Column(Float)
    usertime          = Column(Float)
    elapsed           = Column(String)
    cpu               = Column(Float)
    max_resident_time = Column(Float)

    def __repr__(self):
        return "<CustomToolTypeMixin(systime='%d', usertime='%d', elapsed='%s', cpu='%d', max_resident_time='%d'>" % \
                (self.systime, self.usertime, self.elapsed, self.cpu, self.max_resident_time)

class CustomToolMd5TypeMixin(object):
    ''' Gather timing metrics with input/output uuids and md5 '''
    id                = Column(Integer, primary_key=True)
    case_id           = Column(String)
    vcf_id            = Column(String)
    src_vcf_id        = Column(String)
    tool              = Column(String)
    files             = Column(ARRAY(String))
    systime           = Column(Float)
    usertime          = Column(Float)
    elapsed           = Column(String)
    cpu               = Column(Float)
    max_resident_time = Column(Float)
    md5               = Column(String)

    def __repr__(self):
        return "<CustomToolTypeMd5Mixin(systime='%d', usertime='%d', elapsed='%s', cpu='%d', max_resident_time='%d'>" % \
                (self.systime, self.usertime, self.elapsed, self.cpu, self.max_resident_time)
