
# LEOrbit Python library

This package is aimed to track a satellite in low orbit around Earth. 
It has been created in order to help small space centers (from universities for example) to do calculations for their cubesats. 
But you can use it in any way you want!

### Quick links
- Documentation & reference: ...

## Quick start

### Installing

#### From PyPi

LEOrbit will be available on `PyPi` very soon.

#### Locally (for development purposes)

```sh
git clone -b master https://github.com/grdleo/leorbit.git <leorbit-dir>
cd <leorbit-dir>
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

### Main depedencies
Uses `numpy` internally for computation purposes.

## Quick tour

```python
>>> from leorbit.api import get_satellite, get_passes, Timestamp, Quantity, GPS, VisibleFromEarthLocationEvent, TimeInterval

# Let's compute the position of a LEO satellite!

>>> iss = get_satellite(25544) # 25544: ISS NORAD Cat ID
>>> now = Timestamp.now()
>>> c = iss.coordinates(now) # We compute the ISS position at given time
>>> c.gps()
'<GPS:  003° 59′ 29″E,  042° 29′ 24″S>'

# Now that we have computed its position for the given time, we can project
# it into any frame we want!

# Coordinates in GCRF!
>>> c.gcrf().human_repr("km")
'<Vector3 x=-6381.3991361611925 y=146.7020812274849 z=-2333.735440753873 [km]>'

# Coordinates in ITRF!
>>> c.itrf().human_repr("km")
'<Vector3 x=-2639.11482315943 y=5811.957448727188 z=-2333.735440753873 [km]>'

# Horizontal coordinates in any Earth local frame!
# For example, let's try in Paris.
>>> gps_paris = GPS(
    longitude=2.333333 * Quantity.degree, 
    latitude=48.866667 * Quantity.degree, 
    altitude=0 * Quantity.meter
)
>>> c.horizontal(gps_paris.earth_local_frame)
'<Horizontal: Azimuth:  087° 21′ 36″, Altitude: - 058° 36′ 26″>'

# Ugh, negative altitude, it means it is not visible currently...
# Want to get the ISS passes for the next 7 days?

>>> timeline = TimeInterval(
    start=now,
    stop=now + 7 * Quantity.day,
    dt=5 * Quantity.second
)
>>> passes = get_passes(iss, timeline, gps_paris)

...

```

## Examples & documentation

LEOrbit is available with a complete walkthrough the capabilities of the library.
Check the documentation with the available examples.

# What's next?
- SDP4 implementation (see https://github.com/Bill-Gray/sat_code)
- Orientation tracking for satellites
- Tools to generate a list of commands to send to satellite

### Links
- [Documentation & reference](https://leorbit.readthedocs.org)
- [Developer's page (Léo G.)](https://leog.dev)