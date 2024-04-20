"""Python library `leorbit` makes LEO satellites tracking easy: 
GP data gathering, orbit propagation, coordinates conversions, events computation...

## Reference
This the reference of the library. A complete documentation is in progress.

## Tutorial
[Here](https://github.com/grdleo/leorbit/tree/master/examples) are a full Jupyter notebook walkthrough, as long as a few examples scripts.

## Links
- [Github repository](https://github.com/grdleo/leorbit)
- [Developer's page (Léo G.)](https://leog.dev)
"""

__pdoc__ = dict()

from leorbit.math import PintQuantity, Q_
from leorbit.math import Vec3, Coordinates, AbsoluteFrame, RelativeFrame, GPSCoordinates, HorizontalCoordinates
from leorbit.orbit import Satellite, Star, Sun, Moon, OrbitalElements, Propagator, NoPropagator, SGP4
from leorbit.simulation import AstroEvent, InNorthenAuroraArea, VisibleFromLocation, NightTime, BrightnessBelow
from leorbit.time import Time, Timeline, get_intersections_timelines

def get_sat(catnr: int, propagator: Propagator = SGP4, log=False) -> Satellite:
    """Fetches the latest GP data corresponding to given `catnr`, 
    and returns the corresponding `Satellite` object with given propagator (default: SGP4)"""
    oe: OrbitalElements = OrbitalElements.from_celestrak_norad_cat_id(catnr, log=log)
    return Satellite(propagator(oe), oe.name, oe.norad_cat_id)

# DOC: hide
__pdoc__["leorbit.ext"] = False