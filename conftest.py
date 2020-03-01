import pytest
from Bio import AlignIO
import tests.files.input_dicts as input

file_folder = "./tests/files/"


@pytest.fixture()
def ABC_3_seqs_1_region():
    graphbound = gb.GraphBound(ABC_3_seqs_1_region)
    graphbound.create_alignment()

@pytest.fixture()
def ABC_simple_3_region():
    rgb = gb.RegionGraphBound(input.rd, input.od)
    return rgb