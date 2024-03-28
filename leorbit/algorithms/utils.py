"""Special functions with special purposes. Should not be useful for the average user.
"""

from leorbit.math.vector import Vec3
from math import log10, sin, cos
from leorbit.math import HALF_REV, Q_

RADIIE_AA = Q_(6_378_137, "m")**2
RADIIE_A4 = RADIIE_AA**2
RADIIE_BB = Q_(6_356_752, "m")**2
RADIIE_B4 = RADIIE_BB**2

def geocentric_radius_earth(latitude: Q_ | float) -> Q_:
    """Returns the mean radius of Earth at given latitude.
    Earth is considered as a spheroid. 
    
    Algorithm from: https://en.wikipedia.org/wiki/Earth_radius#Geocentric_radius"""
    if isinstance(latitude, Q_):
        assert latitude.check("rad")
    cc = cos(latitude)
    ss = sin(latitude)
    rr = (RADIIE_A4 * cc + RADIIE_B4 * ss) / (RADIIE_AA * cc + RADIIE_BB * ss)
    return rr**.5

def apparent_magnitude(sun: Vec3, sat: Vec3, observer: Vec3, std_mag: float) -> float:
    """Returns the apparent magnitude of a satellite from its standard magnitude,
    to a given observer.

    All positions must be given in the same orthogonal frame.
    
    Algorithm from: https://astronomy.stackexchange.com/questions/28744/calculating-the-apparent-magnitude-of-a-satellite"""
    dir_obs = observer - sat
    dir_sun = sun - sat
    dist_sat: Q_ = abs(dir_obs)
    phase = dir_sun.angle(dir_obs)

    phi_term = sin(phase) + (HALF_REV - phase).m_as("rad") * cos(phase)

    return std_mag + 5 * log10(dist_sat.m_as("megameter")) - 2.5 * log10(phi_term)

def humanize_duration(t: Q_) -> str:
    """Make a `pint.Quantity` duration human-readable.
    
    Example:
    --------
    
    ```python
    humanize_duration(Q_("723459.23min"))
    >>> '1 year 4 month 15 day 9 hour 39 min 13 s'
    ```
    """
    try:
        assert t.check("s")
    except:
        ValueError(f"{t} is not a `pint.Quantity` or does not have the 'duration' dimension.")

    stages = ["year", "month", "day", "hour", "min", "s"]
    components = {s: 0 for s in stages}

    for s in stages:
        stage_d = Q_(f"1{s}")
        if t < stage_d:
            continue
        c = int(t.m_as(s))
        components[s] = c
        t -= stage_d * c
            
    
    return " ".join(f"{v} {k}" for k, v in components.items() if v > 0)
    


