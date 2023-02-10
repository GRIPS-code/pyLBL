import pytest

from pyLBL import TipsWebApi
from pyLBL.tips import TotalPartitionFunction


def test_tips():
    t, data = TipsWebApi().download("H2O")
    tips = TotalPartitionFunction("H2O", t, data)
    value = tips.total_partition_function(279.54, 1)
    assert value == pytest.approx(160.2790023803711)
