from pyLBL import continua, cross_sections, models, molecular_lines


def test_plugins():
    assert set(molecular_lines.keys()) == set(["arts", "pyLBL"])
    assert set(continua.keys()) == set(["mt_ckd"])
    assert set(cross_sections.keys()) == set(["arts_crossfit"])
    assert set(models) == set(["arts", "arts_crossfit", "mt_ckd", "pyLBL"])
