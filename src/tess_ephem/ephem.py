"""Defines the main interface, i.e. the `ephem` function."""

from functools import lru_cache
from typing import Optional

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
from astroquery.jplhorizons import Horizons
from pandas import DataFrame, concat
from scipy.interpolate import CubicSpline
from tesswcs import locate, pointings

from . import log
from .angle import create_angle_interpolator


class TessEphem:
    def __init__(
        self,
        target: str,
        start: Time = "2018-05-01",
        stop: Time = "2021-01-01",
        step: str = "12H",
        location: str = "@TESS",
        id_type: str = "smallbody",
        orbital_elements: bool = False,
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
        orbital_elements : bool
            If True, returns osculating orbital elements for target (eccentricity, inclination, perihelion distance).
        """
        log.info(f"Started querying JPL/Horizons for ephemeris (id='{target}')")
        eph = _get_horizons_ephem(
            id=target,
            id_type=id_type,
            start=start,
            stop=stop,
            step=step,
            location=location,
        )
        self._raf = create_angle_interpolator(
            eph["datetime_jd"], eph["RA"], enforce_positive=True
        )
        self._decf = create_angle_interpolator(eph["datetime_jd"], eph["DEC"])
        if "V" in eph.columns:
            mag = eph["V"]
        else:
            mag = eph["Tmag"]  # total comet magnitude
        self._vf = CubicSpline(eph["datetime_jd"], mag)
        # Absolute magnitude
        if "H" in eph.columns:
            self._Hf = CubicSpline(eph["datetime_jd"], eph["H"])
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

        if orbital_elements:
            log.info(
                f"Started querying JPL/Horizons for orbital elements (id='{target}')"
            )
            self.elements = _get_horizons_elements(
                id=target,
                id_type=id_type,
                start=start,
                stop=stop,
                step=step,
            )

    def predict_sky(self, time: Time) -> DataFrame:
        ra = self._raf(time.jd)
        dec = self._decf(time.jd)
        v = self._vf(time.jd)
        # If target does not have H magnitude, set to nan.
        if hasattr(self, "_Hf"):
            H = self._Hf(time.jd)
        else:
            H = np.full(len(time.jd), np.nan)
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
                "hmag": H,
                "sun_distance": r,
                "obs_distance": delta,
                "phase_angle": phi,
            }
        )

    def predict(self, time: Time) -> DataFrame:
        """
        Predicts position of target at times.

        Parameters
        ----------
        time : Time
            Times for which to compute the ephemeris.

        Returns
        -------
        ephemeris : DataFrame
            One row for each time stamp that matched a TESS observation.
        orbital_elements_dict : dict
            Average perihelion distance [AU], eccentricity and orbital inclination [deg] of the target during the queried time.
            This is only returned if orbital_elements = True.
        """

        if not isinstance(time, Time):
            time = Time(time)
        if time.isscalar:
            time = time.reshape((1,))

        sky = self.predict_sky(time)
        crd = SkyCoord(sky.ra, sky.dec, unit="deg")

        # Determine TESS sector corresponding to each time
        sectors = []
        for t in time.jd:
            sectors.append(
                pointings["Sector"].value[
                    np.logical_and(t > pointings["Start"], t < pointings["End"])
                ]
            )

        log.info("Started matching the ephemeris to TESS observations")
        # Get sector, camera, ccd, col, row
        df = DataFrame()
        for s in np.unique(np.hstack(sectors)):
            sector_mask = [s in sublist for sublist in sectors]
            result = locate.get_pixel_locations(crd[sector_mask], sector=s).to_pandas()
            result["time"] = time[sector_mask][result["Target Index"]]
            df = concat(
                [df, result[["time", "Sector", "Camera", "CCD", "Column", "Row"]]]
            )

        if len(df) == 0:
            log.warning("Warning: Target not observed by TESS at defined times.")
        else:
            # Reorder df in ascending time
            df = df.sort_values(by="time", ascending=True)
            # Make column names lowercase in df
            df.columns = [x.lower() for x in df.columns]
            # Make sector, camera, ccd integers
            df = df.astype({"sector": int, "camera": int, "ccd": int})
            df = df.merge(sky, on="time", how="inner")
            df = df.set_index("time")

        if hasattr(self, "elements"):
            return df, {
                "perihelion_distance": np.nanmean(self.elements["q"]),
                "eccentricity": np.nanmean(self.elements["e"]),
                "orbital_inclination": np.nanmean(self.elements["incl"]),
            }

        else:
            return df

    @staticmethod
    def from_sector(
        target: str,
        sector: Optional[int] = None,
        id_type: str = "smallbody",
        step: str = "12H",
        return_time: bool = False,
        orbital_elements: bool = False,
    ):
        """
        Initialises TessEphem object from sector number.

        Parameters
        ----------
        target : str
            Horizons target ID.
        sector : int or None
            Sector number.
        id_type : str
            JPL/Horizons target identifier type.
            One of "smallbody", "majorbody", "designation", "name", "asteroid_name",
            "comet_name", or "designation".
        step : str
            Resolution at which ephemeris data will be obtained from JPL Horizons.
        return_time: bool
            If True then return start and stop time of sector.
        orbital_elements : bool
            If True, returns osculating orbital elements for target (eccentricity, inclination, perihelion distance).

        Returns
        -------
        TessEphem
            Initialiased TessEphem object.
        start : Time
            If return_time, sector start time.
        stop : Time
            If return_time, sector end time.
        """

        # If no sector passed, use all sectors in pointings file.
        if sector is None:
            # Define start and stop using pointings file from tesswcs
            start = Time(pointings["Start"][0], format="jd")
            stop = Time(pointings["End"][-1], format="jd")

        else:
            # Sector must be in pointings
            if sector not in pointings["Sector"]:
                raise KeyError(
                    "Sector must be in range {0}-{1}".format(
                        pointings["Sector"][0], pointings["Sector"][-1]
                    )
                )

            # Define start and stop of sector using pointings file from tesswcs
            start = Time(
                pointings[pointings["Sector"] == sector]["Start"][0], format="jd"
            )
            stop = Time(pointings[pointings["Sector"] == sector]["End"][0], format="jd")

        # Buffer of one day on start, stop for future interpolation of ephemeris.
        start_buffer = start - TimeDelta(1, format="jd")
        stop_buffer = stop + TimeDelta(1, format="jd")

        # Return TessEphem object (and start/stop, if return_time=True)
        if return_time:
            return (
                TessEphem(
                    target,
                    start=start_buffer,
                    stop=stop_buffer,
                    step=step,
                    id_type=id_type,
                    orbital_elements=orbital_elements,
                ),
                start,
                stop,
            )
        else:
            return TessEphem(
                target,
                start=start_buffer,
                stop=stop_buffer,
                step=step,
                id_type=id_type,
                orbital_elements=orbital_elements,
            )


@lru_cache()
def _get_horizons_ephem(
    id,
    start: Time,
    stop: Time,
    step: str = "12H",
    id_type: str = "smallbody",
    location: str = "@TESS",
    quantities: str = "1,3,9,19,20,43",
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


@lru_cache()
def _get_horizons_elements(
    id,
    start: Time,
    stop: Time,
    step: str = "12H",
    id_type: str = "smallbody",
):
    """Returns JPL Horizons orbital elements.

    This is simple cached wrapper around astroquery's Horizons.elements.

    By not specifying a location, the heliocenter is used by default.
    """
    epochs = {"start": start.iso, "stop": stop.iso, "step": step}
    log.debug(f"Horizons query parameters:\n\tid={id}\n\tepochs={epochs}")
    t = Horizons(id=id, id_type=id_type, epochs=epochs)
    result = t.elements()
    log.debug(f"Received {len(result)} elements results")
    return result


def ephem(
    target: str,
    time: Time = None,
    sector: Optional[int] = None,
    time_step: float = 1.0,
    id_type: str = "smallbody",
    interpolation_step: str = "12H",
    orbital_elements: bool = False,
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
    time_step : float
        Resolution of time grid if ``sector`` is passed, in days. Will be ignored if ``time`` is passed.
    id_type : str
        JPL/Horizons target identifier type.
        One of "smallbody", "majorbody", "designation", "name", "asteroid_name",
        "comet_name", or "designation".
    interpolation_step : str
        Resolution at which ephemeris data will be obtained from JPL Horizons.
    orbital_elements : bool
        If True, returns osculating orbital elements for target (eccentricity, inclination, perihelion distance).

    Returns
    -------
    ephemeris : DataFrame
        One row for each time stamp that matched a TESS observation.
    orbital_elements_dict : dict
        Average perihelion distance [AU], eccentricity and orbital inclination [deg] of the target during the queried time.
        This is only returned if orbital_elements = True.
    """
    if time is None:
        te, start, stop = TessEphem.from_sector(
            target,
            sector,
            step=interpolation_step,
            id_type=id_type,
            return_time=True,
            orbital_elements=orbital_elements,
        )

        # Define time array using start, stop and time_step.
        days = np.ceil((stop - start).sec / (60 * 60 * 24))
        time = start + TimeDelta(np.arange(0, days + time_step, time_step), format="jd")

    else:
        if not isinstance(time, Time):
            time = Time(time)
        if time.isscalar:
            time = time.reshape((1,))
        # Buffer of seven days on start, stop for future interpolation of ephemeris.
        # (Larger buffer than when defining start, stop from_sector() because, in this case, time might be a single value.)
        start = time[0] - TimeDelta(7, format="jd")
        stop = time[-1] + TimeDelta(7, format="jd")
        te = TessEphem(
            target,
            start=start,
            stop=stop,
            step=interpolation_step,
            id_type=id_type,
            orbital_elements=orbital_elements,
        )
    return te.predict(time=time)
