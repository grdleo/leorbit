
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
Uses `numpy` internally for computation purposes. The consequence is that any computation is absurdly fast.

## Quick tour

Get the coordinates of the satellite of your choice, at the time of your choice
```python
from leorbit.api import get_satellite, get_passes, Timestamp, Qty, GPS, VisibleFromEarthLocationEvent, TimeInterval, N, E

iss = get_satellite(25544) # 25544: ISS NORAD Cat ID
now = Timestamp.now()
coords = iss.coordinates(now) # ISS Coordinates at given time, frame agnostic!

coords.gps()
'<GPS:  003° 59′ 29″E,  042° 29′ 24″S>'
```

Need the coordinates in a specific frame? Sure, as easy as this:
```python
coords.gcrf().human_repr("km") # Coordinates in GCRF!
'<Vector3 x=-6381.3 y=146.7 z=-2333.7 [km]>'

coords.itrf().human_repr("km") # Coordinates in ITRF!
'<Vector3 x=-2639.1 y=5811.9 z=-2333.7 [km]>'
```

Even local coordinates are astonishingly easy to convert to.
```python
gps_paris = 2.333333 * E + 48.866667 * N

coords.horizontal(gps_paris.earth_local_frame)
'<Horizontal: Azimuth:  087° 21′ 36″, Altitude: - 058° 36′ 26″>'
```

You need to know when your satellite passes above your location?
Let's compute the passes for the next 7 days.
```python
timeline = TimeInterval(
    start=now,
    stop=now + 7 * Qty.day,
    dt=5 * Qty.second
)

first_pass, *others = get_passes(iss, timeline, gps_paris)
'<TimeInterval from: \'2026-04-28 at 00:19:19\' to: \'2026-04-28 at 00:28:54\' dt: 5s>'
```

And for the interval of your choice, generate an export of the trajectory:
```python
iss.trajectory(first_pass).horizontal(gps_paris.earth_local_frame).to_csv(first_pass)
"""
timestamp,azimuth,elevation,range
2026-04-28T00:19:21.1,-159.7,0.3,6790.4
2026-04-28T00:19:26.1,-160.1,0.5,6790.4
2026-04-28T00:19:31.1,-160.6,0.8,6790.4
2026-04-28T00:19:36.1,-161.0,1.1,6790.3
...
"""
```

## Examples & documentation

The above examples are one of the few features the library offers.
Make sure to check the documentation and reference to know everything it can do!

LEOrbit is available with a complete walkthrough the capabilities of the library.

## What's next?

The following are features I would like to implement, one day, when I have the time to do so.

- Ephemerids for planets/stars
- More "ready to use out-of-the-box" events, like day/night time, transits, etc.
- SDP4 implementation (see https://github.com/Bill-Gray/sat_code)

## Please participate

You like this library and feel like it is missing a feature? By all means, feel free to open a PR!

## Links
- [Documentation & reference](https://leorbit.readthedocs.org)
- [Developer's page (Léo G.)](https://leog.dev)