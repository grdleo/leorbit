from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from leorbit.propagator import Propagator
    from leorbit.sky_object import Satellite

def get_satellite(
    norad_cat_id: int, 
    propagator: type["Propagator"] | None = None,
    log: bool = True
) -> "Satellite":
    """Returns a `Satellite` object with latest GP data from given `norad_cat_id`, using Celestrak.org database.
    Default propagator used is `SGP4`."""
    
    from leorbit.ext import get_celestrak_gpdata
    from leorbit.propagator import SGP4
    from leorbit.sky_object import Satellite

    if propagator is None:
        propagator = SGP4

    gp_data = get_celestrak_gpdata(norad_cat_id, log)

    return Satellite(
        gp_data.name,
        gp_data.to_orbital_elements(),
        propagator
    )