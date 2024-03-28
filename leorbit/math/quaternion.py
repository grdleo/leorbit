"""FIXME: work in progress"""

from leorbit.math.vector import Vec3
from leorbit.math import Q_

from typing import Any, TypeVar
Q = TypeVar("Q", bound="Quaternion")


class Quaternion:
    def __init__(self, w: float, x: float, y: float, z: float):
        """Creates a quaternion q = w + xi + yj + zk"""
        self._s = w
        self._v = Vec3(x, y, z, _bypass = True)
    
    @staticmethod
    def init_vector_scalar(s: float, v: Vec3):
        return Quaternion(s, float(v._x), float(v._y), float(v._z))
    
    @staticmethod
    def _crash_op_supported(gotten: Any, supported_types: list) -> Exception:
        """Raises an error if operation not supported"""
        supported = [gotten is t for t in supported_types]
        if not any(supported):
            raise ValueError(
                f"Operation not supported with `Vec3` and `{gotten}`."
                f"expected `{supported_types}`"
            )
    
    @property
    def scalar(self: Q) -> float:
        return self._s
    
    @property
    def vector(self: Q) -> Vec3:
        return self._v
    
    def __add__(self: Q, other: Q) -> Q:
        self._crash_op_supported(other, [Q])
        return Quaternion.init_vector_scalar(
            self._s + other._s,
            self._v + other._v
        )
    
    def __neg__(self: Q) -> Q:
        return Quaternion.init_vector_scalar(-self._s, -self._v)
    
    def __sub__(self: Q, other: Q) -> Q:
        return self + -other
    
    def __mul__(self: Q, other: float) -> Q:
        """Scalar multiplication"""
        Quaternion._crash_op_supported(other, [float])
        return Quaternion.init_vector_scalar(self._s * other, self._v * other)
    
    def quat_prod(self: Q, other: Q) -> Q:
        """Quaternion product (self x other)"""
        Quaternion._crash_op_supported(other, [Q])
        s = self._s * other._s - self._v @ other._v
        v = self._v * other._s + other._v * self._s + self._v & other._v
        return Quaternion.init_vector_scalar(s, v)
    
    # q @ r
    def __matmul__(self: Q, other: Q) -> Q:
        return self.quat_prod(other)
    
    def inv(self: Q) -> Q:
        den = self._s**2 + Vec3.dot(self._v, self._v)
        return Quaternion.init_vector_scalar(self._s / den, -self._v / den)
        
        
    
            
