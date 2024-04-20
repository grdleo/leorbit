
# LEOrbit Python library

This package is aimed to track a satellite in low orbit around Earth. It has been created in order to help small space centers (from universities for example) to do calculations for their cubesats. But you can use it in any way you want!

### Reference & documentation
- Reference: https://grdleo.github.io/leorbit

## Quick start

### Installing

LEOrbit will be available on `PyPi` very soon.

The build is already available on [`test.pypi.org`](https://test.pypi.org/project/leorbit/)

### Main depedencies
Uses `pint` to handle physical quantities with units. `requests` and `beautifulsoup4` to fetch data on Celestrak and others.

No NumPy or any other heavy library, `leorbit` wants to be and remain pure Python.

## Capabilities

* **Position prediction** at any epoch of any satellite in LEO (low-earth orbit) with a NORAD catalog ID, or manually fetched GP data, or even state vectors.

```python
>>> from leorbit import get_sat, Time, GPSCoordinates, Timeline, VisibleFromLocation, Q_

# Creates a `Satellite` object from latest GP data, ready to be propagated using SGP4
>>> iss = get_sat(25544) 
<Satellite name='ISS (ZARYA)' t₀='2024-03-19T08:52:21.888192+00:00' [SGP4]>

>>> t = Time.now()
>>> c = iss.coordinates(t) # Position at epoch `t` stored here
>>> c.to_gps()
GPSCoordinates(longitude_deg=-90.5461239, latitude_deg=-43.8902309)
```

* **Astonishingly easy coordinates conversions** from/to any frame, like ITRF, GCRF, or any relative frame created by user.

```python
# Coordinates in GCRF 
>>> c.to_gcrf()
Vec3(x=4813.582 km, y=-918.576 km, z=-4714.194 km)

# Coordinates in ITRF 
>>> c.to_itrf()
Vec3(x=-46.708 km, y=-4900.222 km, z=-4714.194 km)

# Horizontal coordinates (location=Paris)
>>> gps_paris = GPSCoordinates(2.333333, 48.866667, 0., "Paris (FR)")
>>> c.to_horizontal_coordinates(gps_paris)
HorizontalCoordinates(azimuth_deg=-120.7834, altitude_deg=-60.5245)
```

* **Computing intervals where satellite meets conditions** like
    * Being visible from a certain location on Earth
    * Having a certain brightness to an observer
    * Satellite enlighten, Moon visible, ...
    * ... or any other event created by user
    * ... or any combinations of previous event using `intersections` and `unions`

```python
>>> now = Time.now()
# Full timeline of computation
>>> simu_tl = Timeline(now, now + Q_("7 days"), dt=Q_("5s"))

# List of time intervals when ISS is visible in Paris
>>> visi_intervals = VisibleFromLocation(iss, gps_paris).compute_intervals(simu_tl) 
>>> [interval.human for interval in visi_intervals]
[
    "Timeline: from '2024-03-22T22:21:04' to '2024-03-22T22:26:30', duration=5 min 25 s",
    "Timeline: from '2024-03-24T03:06:25' to '2024-03-24T03:11:11', duration=4 min 45 s"
]
```

## Examples & documentation

LEOrbit is available with a complete walkthrough the capabilities of the library.
It also have many examples all available [here](https://github.com/grdleo/leorbit/tree/master/examples).

# What's next?
- SDP4 implementation (see https://github.com/Bill-Gray/sat_code)
- Joblib implementation for parallel computations https://joblib.readthedocs.io/en/latest/
- Orientation tracking for satellites
- Tools to generate a list of commands to send to satellite

### Links
- [Reference & documentation](https://grdleo.github.io/leorbit)
- [Developer's page (Léo G.)](https://leog.dev)