tess-ephem
==========

**Where are Solar System objects located in TESS FFI data?***

`tess-ephem` is a user-friendly package which enables users to compute the positions of Solar System objects -- asteroids, comets, and planets --
in the data archive of NASA's TESS Space Telescope.

Installation
------------

.. code-block:: bash

    python -m pip install tess-ephem


Example use
-----------

tess-ephem allows you to search the entire archive of TESS FFI's for the presence
of a known minor planet, and obtain the result as a Pandas DataFrame.
For example:

.. code-block:: python

    >>> from tess_ephem import ephem
    >>> ephem("Sedna")
                             sector  camera  ccd       column          row
    time
    2018-11-16 00:00:00.000       5       1    4  1543.312296  1102.821559
    2018-11-17 00:00:00.000       5       1    4  1545.160910  1102.880825
    2018-11-18 00:00:00.000       5       1    4  1547.011351  1102.934375
    ...
    2018-12-09 00:00:00.000       5       1    4  1584.585407  1102.239292
    2018-12-10 00:00:00.000       5       1    4  1586.245261  1102.132304
    2018-12-11 00:00:00.000       5       1    4  1587.906380  1102.012091


You can also obtain the ephemeris for one or more specific times
by passing the `time` parameter:

.. code-block:: python

    >>> ephem("Sedna", time="2018-11-21 17:35:00")
                             sector  camera  ccd       column          row
    time
    2018-11-21 17:35:00.000       5       1    4  1553.887838  1103.048431


Additional physical parameters can be obtained by passing the `verbose=True` parameter:

.. code-block:: python

    >>> ephem("Sedna", time="2018-11-21 17:35:00", verbose=True)
                             sector  camera  ccd       column          row  pixels_per_hour        ra      dec    vmag  sun_distance  obs_distance  phase_angle
    time
    2018-11-21 17:35:00.000       5       1    4  1553.887838  1103.048431         0.074054  57.05786  7.63721  20.612     84.942885     83.975689       0.1419


Finally, using the companion `tess-locator <https://github.com/SSDataLab/tess-locator>`_
package, you can convert the TESS pixel coordinates directly to FFI filenames and urls:

.. code-block:: python

    >>> from tess_locator import TessCoordList
    >>> df = ephem("Sedna", time=["2018-11-21 17:35:00", "2018-11-22 17:35:00"])
    >>> TessCoordList.from_pandas(df).get_images()
    List of 2 images
    ↳[TessImage(filename='tess2018325165939-s0005-1-4-0125-s_ffic.fits', begin='2018-11-21 17:07:46', end='2018-11-21 17:37:46')
      TessImage(filename='tess2018326165939-s0005-1-4-0125-s_ffic.fits', begin='2018-11-22 17:07:46', end='2018-11-22 17:37:46')]
