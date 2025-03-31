from copy import copy
from typing import Self
from mathematics.mat33 import Mat33, MAT_I
from mathematics.vec3 import Vec3, VEC_0, POS_UNIT
from functools import reduce

class Transform:
    """A transform that will change a `Vec3`
    """
    def apply(self, v: Vec3) -> Vec3:
        raise NotImplementedError()
    
    def unapply(self, v: Vec3) -> Vec3:
        raise NotImplementedError()
    
    def reverse(self) -> Self:
        tr = copy(self)
        tr.apply = self.unapply
        tr.unapply = self.apply
        return tr

class Identity(Transform):
    def __init__(self):
        pass

    def apply(self, v: Vec3) -> Vec3:
        return v
    
    def unapply(self, v: Vec3) -> Vec3:
        return v
    
class LinearTransform(Transform):
    def __init__(self, matrix: Mat33):
        self.matrix = matrix

    def apply(self, v: Vec3) -> Vec3:
        return self.matrix.prod_vec(v)
    
    def unapply(self, v: Vec3) -> Vec3:
        return self.matrix.inverse.prod_vec(v)

class AffineTransform(Transform):
    def __init__(self, linear_transform: Mat33 = MAT_I, translation: Vec3 = VEC_0):
        self.linear_transform = linear_transform
        self.translation = translation
    
    def apply(self, v: Vec3) -> Vec3:
        """Applies the matrix to the given vector, then the translation (if is a position)."""
        vv = self.linear_transform.prod_vec(v)
        if v.unit.is_compatible_with(POS_UNIT):
            vv += self.translation
        return vv
    
    def unapply(self, v: Vec3) -> Vec3:
        """Returns the vector that results the given vector when this transformation is applied."""
        vv = (v - self.translation) if v.unit.is_compatible_with(POS_UNIT) else v
        return self.linear_transform.inverse.prod_vec(vv)
    
class ChainTransform(Transform):
    def __init__(self, *transforms: Transform):
        """First applied first"""
        self.transforms: list[Transform] = [t for t in transforms if not isinstance(t, Identity)]
    
    def apply(self, v: Vec3) -> Vec3:
        if not self.transforms:
            return v
        return reduce(lambda _v, transform: transform.apply(_v), self.transforms, v)
    
    def unapply(self, v: Vec3) -> Vec3:
        if not self.transforms:
            return v
        return reduce(lambda _v, transform: transform.unapply(_v), reversed(self.transforms), v)

    def compose(self, *transforms: AffineTransform) -> Self:
        """Returns a new transformation composed with given single transform.
        Given transform will be performed last."""
        return ChainTransform(*self.transforms, *transforms)
    
    def compose_chain(self, chain: Self) -> Self:
        """Returns a new transformation composed with given transform.
        Given transform will be performed last."""
        return ChainTransform(*self.transforms, *chain.transforms)

### SPECIAL TRANSFORMS ###

class RotationZ(Transform):
    """A simple transform that represent a rotation around the `z` axis."""
    def __init__(self, angle: float):
        """`angle [rad]`"""
        self.angle = angle

    def apply(self, v: Vec3) -> Vec3:
        return v.rotate_zaxis(self.angle)
    
    def unapply(self, v: Vec3) -> Vec3:
        return v.rotate_zaxis(-self.angle)