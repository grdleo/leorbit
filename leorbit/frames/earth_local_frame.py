from copy import copy
from math import cos
from frames.absolute_frame import AbsoluteFrame
from frames.relative_frame import RelativeFrame
from mathematics.mat33 import Mat33
from mathematics.transformation import AffineTransform, LinearTransform
from mathematics.vec3 import Vec3, ZAXIS, UREG
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from coordinates.coordinates import Coordinates

class EarthLocalFrame(RelativeFrame):
    """Coordinates frame relative to a given location on Earth.
    
        - Origin: Given location
        - `z:` towards zenith (aka. towards the sky, perpendicular to ground)
        - `y:` towards "East"
        - `x × y = -z`
    """
    location: "Coordinates"

    def __init__(self, location: "Coordinates"):
        itrf: Vec3 = location.get_pos(AbsoluteFrame.ITRF)
        
        z = itrf.normalized() # towards zenith
        north = copy(ZAXIS)
        ang = Vec3.angle(z, north)
        x: Vec3 # towards "true north"

        if ang % (180 * UREG.degrees) == 0: # FIXME
            raise ValueError("Cannot create `EarthLocalFrame` in Earth's poles!")
        elif ang == 90 * UREG.degrees: # FIXME
            x = copy(north)
        else:
            x = (north / cos(ang) - z).normalized()
            if ang > 90 * UREG.degrees:
                x *= -1
        
        y = x.cross(z) # towards "east"

        mat = Mat33.from_col_vectors(x, y, z)
        transform = AffineTransform(mat, itrf)

        super().__init__(AbsoluteFrame.ITRF, transform)

        self.location = location
    
    def __repr__(self) -> str:
        gps = self.location.gps()
        return f"<EarthLocalFrame at GPS location {gps.dms}>"