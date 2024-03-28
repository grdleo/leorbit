"""Special functions with special purposes. Should not be useful for the average user."""

from leorbit.algorithms.utils import apparent_magnitude
from leorbit.math.coordinate import EarthLocalFrame, GPSCoordinates
from leorbit.math.vector import Vec3
from leorbit.orbit import SUN
from leorbit.orbit.objects import Satellite, Sun
from leorbit.time.time import Time

def compute_magnitude(epoch: Time, local: EarthLocalFrame, sat: Satellite, sun: Sun = SUN) -> float:
    """Returns the apparent magnitude of given sat at given location, at given epoch."""
    sun_coord = sun.coordinates(epoch)
    sat_coord = sat.coordinates(epoch)

    sun_vec = sun_coord.pos_in_frame(local)
    sat_vec = sat_coord.pos_in_frame(local)

    return apparent_magnitude(sun_vec, sat_vec, Vec3.zero("m"), sat.std_mag)