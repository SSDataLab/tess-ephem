tess-ephem
==========


**Where are Solar System objects located in TESS FFI data?**

|pypi| |pytest| |black| |flake8| |mypy|

.. |pypi| image:: https://img.shields.io/pypi/v/tess-ephem
                :target: https://pypi.python.org/pypi/tess-ephem
.. |pytest| image:: https://github.com/SSDataLab/tess-ephem/workflows/pytest/badge.svg
.. |black| image:: https://github.com/SSDataLab/tess-ephem/workflows/black/badge.svg
.. |flake8| image:: https://github.com/SSDataLab/tess-ephem/workflows/flake8/badge.svg
.. |mypy| image:: https://github.com/SSDataLab/tess-ephem/workflows/mypy/badge.svg

`tess-ephem` is a user-friendly package which enables users to compute the positions of Solar System objects -- asteroids, comets, and planets --
in the data archive of NASA's TESS Space Telescope.

Installation
------------

.. code-block:: bash

    python -m pip install tess-ephem


Example use
-----------

tess-ephem allows you to search the entire archive of TESS FFI's for the presence
of a known minor planet, and obtain the result as a Pandas DataFrame. The output pixel coordinates (column and row) follow the TESS convention, with (1,1) being the middle of the pixel in the lower left corner of the FFI. For example:

.. code-block:: python

    >>> from tess_ephem import ephem
    >>> ephem("Sedna")
               sector  camera  ccd       column          row
    time                                                    
    2458437.5       5       1    4  1540.328759  1102.742761
    2458438.5       5       1    4  1542.057935  1102.906116
    2458439.5       5       1    4  1543.919678  1102.977150
    2458440.5       5       1    4  1545.806011  1103.011147
    2458441.5       5       1    4  1547.691635  1103.029184
    ...           ...     ...  ...          ...          ...
    2460254.5      71       2    4  1984.472509  1004.531966
    2460255.5      71       2    4  1984.704905  1002.716266
    2460256.5      71       2    4  1984.934016  1000.892089
    2460257.5      71       2    4  1985.160431   999.062904
    2460258.5      71       2    4  1985.376804   997.240991

    [78 rows x 5 columns]


You can also obtain the ephemeris for one or more specific times
by passing the `time` parameter:

.. code-block:: python

    >>> ephem("Sedna", time="2018-11-21 17:35:00")
                             sector  camera  ccd       column          row
    time                                                                  
    2018-11-21 17:35:00.000       5       1    4  1552.813087  1103.033716

    >>> from astropy.time import Time
    >>> ephem("Sedna", time=Time([2458441.5,2460258.5], format='jd'))
               sector  camera  ccd       column          row
    time                                                    
    2458441.5       5       1    4  1547.691635  1103.029184
    2460258.5      71       2    4  1985.376804   997.240991


Additional physical parameters can be obtained by passing the `verbose=True` parameter:

.. code-block:: python

    >>> ephem("Sedna", time="2018-11-21 17:35:00", verbose=True)
                             sector  camera  ccd       column          row  pixels_per_hour        ra      dec    vmag  sun_distance  obs_distance  phase_angle
    time                                                                                                                                                       
    2018-11-21 17:35:00.000       5       1    4  1552.813087  1103.033716         0.074053  57.06362  7.63836  20.812     84.943049     83.975854       0.1419


You can alternatively obtain the ephemeris during a specific sector by passing 
the `sector` parameter:

.. code-block:: python

    >>> ephem("Sedna", sector=70)
               sector  camera  ccd       column          row
    time                                                    
    2460208.5      70       4    2  1965.819900  1827.440280
    2460209.5      70       4    2  1966.122988  1826.880450
    2460210.5      70       4    2  1966.445615  1826.219237
    2460211.5      70       4    2  1966.792833  1825.480366
    2460212.5      70       4    2  1967.156084  1824.685065
    2460213.5      70       4    2  1967.530374  1823.844978
    2460214.5      70       4    2  1967.912846  1822.964230
    2460215.5      70       4    2  1968.300642  1822.046948
    2460216.5      70       4    2  1968.693056  1821.098583
    2460217.5      70       4    2  1969.085076  1820.121939
    2460218.5      70       4    2  1969.477787  1819.122100
    2460219.5      70       4    2  1969.865471  1818.107325
    2460220.5      70       4    2  1970.236706  1817.102989
    2460221.5      70       4    2  1970.537507  1816.171600
    2460222.5      70       4    2  1970.786337  1815.215528
    2460223.5      70       4    2  1971.057940  1814.164426
    2460224.5      70       4    2  1971.352361  1813.044830
    2460225.5      70       4    2  1971.660316  1811.874587
    2460226.5      70       4    2  1971.976449  1810.663652
    2460227.5      70       4    2  1972.300053  1809.417480
    2460228.5      70       4    2  1972.626477  1808.140569
    2460229.5      70       4    2  1972.954292  1806.834984
    2460230.5      70       4    2  1973.282790  1805.506180
    2460231.5      70       4    2  1973.609473  1804.159986
    2460232.5      70       4    2  1973.931842  1802.802230


When passing the `sector` parameter, the `time_step` is by default 1 day. 
This can be changed as follows:

    >>> ephem("Sedna", sector=70, time_step=0.1)
               sector  camera  ccd       column          row
    time                                                    
    2460207.6      70       4    2  1965.495431  1827.937212
    2460207.7      70       4    2  1965.535648  1827.878206
    2460207.8      70       4    2  1965.575019  1827.820108
    2460207.9      70       4    2  1965.613392  1827.763020
    2460208.0      70       4    2  1965.650616  1827.707041
    ...           ...     ...  ...          ...          ...
    2460233.0      70       4    2  1974.086940  1802.125478
    2460233.1      70       4    2  1974.117634  1801.990490
    2460233.2      70       4    2  1974.148118  1801.855903
    2460233.3      70       4    2  1974.178192  1801.721961
    2460233.4      70       4    2  1974.207660  1801.588906

    [259 rows x 5 columns]

