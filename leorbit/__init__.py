from leorbit.propagator import SGP4, Propagator
from leorbit.sky_object import Satellite
from leorbit.coordinates import OrbitalElements
from leorbit.ext import get_celestrak_gpdata

def get_satellite(
    norad_cat_id: int, 
    propagator: type[Propagator] = SGP4,
    log: bool = True
) -> Satellite:
    """Returns a `Satellite` object with latest GP data from given `norad_cat_id`, using Celestrak.org database.
    Default propagator used is `SGP4`."""
    
    gp_data = get_celestrak_gpdata(norad_cat_id, log)

    return Satellite(
        gp_data.name,
        gp_data.to_orbital_elements(),
        propagator
    )