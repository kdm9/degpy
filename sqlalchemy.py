from sqlalchemy import Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
# give all tables a primary key
Base.pk = Column(Integer, primary_key=True)


class TPlot(Base):

    pos = Column(Integer)
    count = Column(Integer)

    def __init__(self):
        super(TPlot, self).__init__()


class Hit(Base):

    read_name = Column(String(50))
    agi_id = Column(Integer)
    pos = Column(Integer)
    cigar = Column(String(50))
    mapq = Column(Integer)
    passed_qc = Column(Integer)
    sample_id = Column(Integer)

    def __init__(self):
        super(Hit, self).__init__()


class Sample(Base):

    def __init__(self):
        super(Sample, self).__init__()


class Target(Base):

    def __init__(self):
        super(Target, self).__init__()


class Histogram(Base):

    def __init__(self):
        super(Histogram, self).__init__()
