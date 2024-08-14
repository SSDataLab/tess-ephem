from tess_ephem import ephem, TessEphem
from astropy.time import Time, TimeDelta
import numpy as np


def test_comet():
    """
    Does `ephem()` work for a comet identifier?
    Horizons provides a column called "Tmag" rather than "V" for a comet's
    total magnitude.  Let's make sure we support this!
    """
    comet_ephem = ephem("90000700", time="2021-08-21", verbose=True)
    assert "vmag" in comet_ephem.columns


def test_from_sector():
    """
    Check from_sector() pulls the correct start/stop time from the pointings.csv file.
    """
    _, start, stop = TessEphem.from_sector("1998 YT6", sector=6, return_time=True)
    assert start.value == 2458463.5
    assert stop.value == 2458490.5


def test_predict():
    """
    Check predict() gives the expected output for an example asteroid.
    """

    # Get pixel location of asteroid "1998 YT6" at time 2458489.8075233004 JD (during sector 6)
    t = Time(2458489.8075233004, format="jd")
    ephem = TessEphem(
        "1998 YT6",
        start=t - TimeDelta(7, format="jd"),
        stop=t + TimeDelta(7, format="jd"),
    )
    pixel_locations = ephem.predict(time=t)

    assert len(pixel_locations) == 1
    assert pixel_locations["sector"][0] == 6
    assert pixel_locations["camera"][0] == 1
    assert pixel_locations["ccd"][0] == 1

    # Expected col, row calculated by:
    # 1. using JPL/Horizons web interface to find RA,Dec of asteroid at time t.
    # 2. using tessrip to get average WCS from sector/camera/ccd
    # 3. converting RA,Dec to col,row with wcs_world2pix()
    assert np.round(pixel_locations["row"][0], 1) == 1107.6
    assert np.round(pixel_locations["column"][0], 1) == 1087.8
