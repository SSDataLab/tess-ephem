"""Defines the main interface, i.e. the `ephem` function."""
from functools import lru_cache
from typing import Optional

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from astroquery.jplhorizons import Horizons
from pandas import DataFrame
from scipy.interpolate import CubicSpline

from tess_locator import locate
from tess_locator.dates import get_sector_dates

from .angle import create_angle_interpolator
from . import log


class TessEphem:
    def __init__(
        self,
        target: str,
        start: Time = "2018-05-01",
        stop: Time = "2021-01-01",
        step: str = "12H",
        location: str = "@TESS",
        id_type: str = "smallbody",
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
            id=target,
            id_type=id_type,
            start=start,
            stop=stop,
            step=step,
            location=location,
        )
        self._raf = create_angle_interpolator(
            eph["datetime_jd"], eph["RA_app"], enforce_positive=True
        )
        self._decf = create_angle_interpolator(eph["datetime_jd"], eph["DEC_app"])
        if "V" in eph.columns:
            mag = eph["V"]
        else:
            mag = eph["Tmag"]  # total comet magnitude
        self._vf = CubicSpline(eph["datetime_jd"], mag)
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

    def predict(
        self, time: Time, aberrate: bool = True, verbose: bool = False
    ) -> DataFrame:
        sky = self.predict_sky(time)
        crd = SkyCoord(sky.ra, sky.dec, unit="deg")
        log.info("Started matching the ephemeris to TESS observations")
        locresult = locate(crd, time=time, aberrate=aberrate)
        df = locresult.to_pandas().merge(sky, on="time", how="inner")
        df = df.set_index("time")
        if not verbose:
            df = df[["sector", "camera", "ccd", "column", "row"]]
        return df


@lru_cache()
def _get_horizons_ephem(
    id,
    start: Time,
    stop: Time,
    step: str = "12H",
    id_type: str = "smallbody",
    location: str = "@TESS",
    quantities: str = "2,3,9,19,20,43",
):
    """Returns JPL Horizons ephemeris.

    This is simple cached wrapper around astroquery's Horizons.ephemerides.
    """
    epochs = {"start": start.iso, "stop": stop.iso, "step": step}
    log.debug(
        f"Horizons query parameters:\n\tid={id}\n\tlocation={location}\n\tepochs={epochs}"
    )
    t = Horizons(id=id, id_type=id_type, location=location, epochs=epochs)
    result = t.ephemerides(quantities=quantities)
    log.debug(f"Received {len(result)} ephemeris results")
    return result


def ephem(
    target: str,
    time: Time = None,
    sector: Optional[int] = None,
    verbose: bool = False,
    id_type: str = "smallbody",
    interpolation_step: str = "12H",
    aberrate: bool = True,
) -> DataFrame:
    """Returns the ephemeris of a Solar System body in the TESS FFI data set.

    Parameters
    ----------
    target : str
        Horizons target ID.
    time : Time
        Times for which to compute the ephemeris.
        By default, one time stamp for every day in the TESS mission will be queried.
    sector : int
        Sector number.  Will be ignored if ``time`` is passed.
    verbose : bool
        Return extra parameters?
    id_type : str
        JPL/Horizons target identifier type.
        One of "smallbody", "majorbody", "designation", "name", "asteroid_name",
        "comet_name", or "designation".
    interpolation_step : str
        Resolution at which ephemeris data will be obtained from JPL Horizons.
    aberrate: bool
        If True then use tess-point's approximate DVA when computing the pixel coordinates.

    Returns
    -------
    ephemeris : DataFrame
        One row for each time stamp that matched a TESS observation.
    """
    if time is None:
        dates = get_sector_dates(sector=sector)
        start = Time(dates.iloc[0].begin[0:10])
        if sector:
            stop = Time(dates.iloc[-1].end[0:10])
        else:
            # Hack: use `-3` because Horizons does not contain TESS ephemeris beyond Jul 2022
            stop = Time(dates.iloc[-3].end[0:10])
        days = np.ceil((stop - start).sec / (60 * 60 * 24))
        time = start + np.arange(-1, days + 1, 1.0)
    else:
        if not isinstance(time, Time):
            time = Time(time)
        if time.isscalar:
            time = time.reshape((1,))
        start = time[0] - 7
        stop = time[-1] + 7
    te = TessEphem(
        target, start=start, stop=stop, step=interpolation_step, id_type=id_type
    )
    return te.predict(time=time, aberrate=aberrate, verbose=verbose)
