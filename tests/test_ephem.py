from tess_ephem import ephem


def test_comet():
    """Does `ephem()` work for a comet identifier?
    Horizons provides a column called "Tmag" rather than "V" for a comet's
    total magnitude.  Let's make sure we support this!
    """
    comet_ephem = ephem("90000700", time="2019-01-01", verbose=True)
    "vmag" in comet_ephem.columns
