"""Defines the main interface, i.e. the `ephem` function."""
from functools import lru_cache

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from astroquery.jplhorizons import Horizons
from pandas import DataFrame
from scipy.interpolate import CubicSpline

import tess_locator as tl
from tess_locator import locate

from .angle import create_angle_interpolator
from . import log


class TessEphem:
    def __init__(
        self,
        target: str,
        start: Time = "2018-05-01",
        stop: Time = "2021-01-01",
        step: str = "7D",
        location: str = "@TESS",
    ):
        """
        Parametrized as start, stop, step because it is the most efficient way to
        get ephemeris from Horizons.

        Parameters
        ----------
        target : str
            JPL/Horizons target ID.
        start : Time
        stop : Time
        step : str
        location : str
        """
        log.info(f"Started querying JPL/Horizons for for ephemeris (id='{target}')")
        eph = _get_horizons_ephem(
            id=target, start=start, stop=stop, step=step, location=location
        )
        self._raf = create_angle_interpolator(
            eph["datetime_jd"], eph["RA"], enforce_positive=True
        )
        self._decf = create_angle_interpolator(eph["datetime_jd"], eph["DEC"])
        self._vf = CubicSpline(eph["datetime_jd"], eph["V"])
        # Sun-target distance
        self._rf = CubicSpline(eph["datetime_jd"], eph["r"])
        # Observer-target distance
        self._deltaf = CubicSpline(eph["datetime_jd"], eph["delta"])
        # Phase angle
        self._phif = CubicSpline(eph["datetime_jd"], eph["alpha_true"])
        # Motion expressed in TESS pixels per hour
        motion = (
            np.hypot(eph["RA_rate"].quantity, eph["DEC_rate"].quantity)
            .to(u.arcsec / u.hour)
            .value
            / 21.0
        )
        self._motionf = CubicSpline(eph["datetime_jd"], motion)
        self.target = target
        self.ephemerides = eph

    def predict_sky(self, time: Time) -> DataFrame:
        ra = self._raf(time.jd)
        dec = self._decf(time.jd)
        v = self._vf(time.jd)
        r = self._rf(time.jd)
        delta = self._deltaf(time.jd)
        phi = self._phif(time.jd)
        motion = self._motionf(time.jd)
        return DataFrame(
            {
                "time": time,
                "pixels_per_hour": motion,
                "ra": ra,
                "dec": dec,
                "vmag": v,
                "sun_distance": r,
                "obs_distance": delta,
                "phase_angle": phi,
            }
        )

    def predict(self, time: Time, verbose: bool = False) -> DataFrame:
        # return self.ephemerides
        sky = self.predict_sky(time)
        crd = SkyCoord(sky.ra, sky.dec, unit="deg")
        log.info("Started matching the ephemeris to TESS observations")
        locresult = locate(crd, time=time)
        df = locresult.to_pandas().merge(sky, on="time", how="inner")
        df = df.set_index("time")
        if not verbose:
            df = df[["sector", "camera", "ccd", "column", "row"]]
        return df


@lru_cache
def _get_horizons_ephem(
    id,
    start: Time,
    stop: Time,
    step="7D",
    location="@TESS",
    quantities="1,3,9,19,20,43",
):
    """Returns JPL Horizons ephemeris.

    This is simple cached wrapper around astroquery's Horizons.ephemerides.
    """
    epochs = {"start": start.iso, "stop": stop.iso, "step": step}
    log.debug(
        f"Horizons query parameters:\n\tid={id}\n\tlocation={location}\n\tepochs={epochs}"
    )
    t = Horizons(id=id, location=location, epochs=epochs)
    result = t.ephemerides(quantities=quantities)
    log.debug(f"Received {len(result)} ephemeris results")
    return result


def ephem(target: str, time: Time = None, verbose: bool = False) -> DataFrame:
    """Returns the ephemeris of a Solar System body in the TESS FFI data set.

    Parameters
    ----------
    target : str
        Horizons target ID.
    time : Time
        Times for which to compute the ephemeris.
        By default, one time stamp for every day in the TESS mission will be queried.

    Returns
    -------
    ephemeris : DataFrame
        One row for each time stamp that matched a TESS observation.
    """
    if time is None:
        dates = tl.wcs_catalog.get_sector_dates()
        start = Time(dates.iloc[0].begin[0:10])
        stop = Time(dates.iloc[-1].end[0:10])
        step = "7D"
        days = np.ceil((stop - start).sec / (60 * 60 * 24))
        time = start + np.arange(-1, days + 1, 1.0)
    else:
        if not isinstance(time, Time):
            time = Time(time)
        if time.isscalar:
            time = time.reshape((1,))
        start = time[0] - 7
        stop = time[-1] + 7
        step = "7D"
    te = TessEphem(target, start=start, stop=stop, step=step)
    return te.predict(time=time, verbose=verbose)
