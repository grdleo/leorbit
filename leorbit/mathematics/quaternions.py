from functools import cached_property

from typing import TYPE_CHECKING, Type

if TYPE_CHECKING:
    from mathematics.vec3 import Vec3
else:
    Vec3 = Type


class Quaternion:
    def __init__(self,
        a: float,
        b: float,
        c: float,
        d: float
    ):
        self._a = a
        self._b = b
        self._c = c
        self._d = d
    
    @staticmethod
    def from_scalar_vector(s: float, v: Vec3) -> "Quaternion":
        return Quaternion(s, v.x, v.y, v.z)

    @property
    def abcd(self) -> tuple[float, float, float, float]:
        return self._a, self._b, self._c, self._d
    
    @property
    def bcd(self) -> tuple[float, float, float]:
        return self._b, self._c, self._d
    
    def quat_mul(self, other: "Quaternion") -> "Quaternion":
        assert isinstance(other, Quaternion)

        sa, sb, sc, sd = self.abcd
        oa, ob, oc, od = other.abcd

        return Quaternion(
            sa*oa - sb*ob - sc*oc - sd*od,
            sb*oa + sa*ob + sc*od - sd*oc,
            sc*oa + sa*oc + sd*ob - sb*od,
            sa*od + sd*oa + sb*oc - sc*ob
        )
    
    def conjugate(self, other: "Quaternion") -> "Quaternion":
        return other @ self @ other.inverse
    
    def __matmul__(self, other: "Quaternion") -> "Quaternion":
        return self.quat_mul(other)

    def __mult__(self, other: float | int) -> "Quaternion":
        assert isinstance(other, float | int)

        return Quaternion(
            self._a * other,
            self._b * other,
            self._c * other,
            self._d * other
        )
    
    def __imul__(self, other: float | int) -> "Quaternion":
        return self.__mul__(other)
    
    def __truediv__(self, other: float | int) -> "Quaternion":
        assert isinstance(other, float | int)

        return self.__mul__(1/other)
    
    @cached_property
    def magnitude_sqr(self) -> float:
        return self._a**2 + self._b**2 + self._c**2 + self._d**2
    
    @cached_property
    def magnitude(self) -> float:
        return self.magnitude_sqr**.5
    
    @cached_property
    def conj(self) -> "Quaternion":
        return Quaternion(self._a, -self._b, -self._c, -self._d)
    
    @cached_property
    def inverse(self) -> "Quaternion":
        return self.conj / self.magnitude_sqr