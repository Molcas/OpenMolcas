import pytest
from FCIDUMP_tools.IntegralClass import *
from GA_mod.GUGA_diag import *
import numpy as np

FILEPATH = "../extra_files/"

diagelem_benchmark = np.load(FILEPATH + "Benzene_diagelem.npy")
csfs = np.load(FILEPATH + "Singlet_6i6_CSFs.npy")

IntegralClass = FCIDUMPReader(FILEPATH + "FCIDUMP_benzene_loc")
DiagElement = DiagElement(6, IntegralClass)


class TestDiagElement:
    @pytest.mark.parametrize("CSF,expected", [
        (csfs[i], diagelem_benchmark[i]) for i in range(len(csfs))
    ])
    def test_diagelements(self, CSF, expected):
        result = DiagElement.calc_diag_elem(CSF)
        assert np.isclose(result, expected)
