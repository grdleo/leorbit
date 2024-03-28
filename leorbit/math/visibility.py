"""FIXME: work in progress"""

from leorbit.math.vector import Vec3
from leorbit.math import Q_
from math import atan, pi, acos

def two_circles_intersection_area(r1: Q_, r2: Q_, d: Q_) -> Q_:
    """Gives the intersection area of two circles of radiuses r1 and r2
    spaced of d"""
    sqr_r1 = r1**2
    sqr_r2 = r2**2
    d1 = 0.5 * (sqr_r1 - sqr_r2 + d**2) / d
    d2 = d - d1

    return (sqr_r1 * acos(d1 / r1) + sqr_r2 * acos(d2 / r2)
        - d1 * (sqr_r1 - d1**2)**0.5 - d2 * (sqr_r2 - d2**2)**0.5)

class Sphere:
    def __init__(self, translated_pos: Vec3, radius: Q_):
        self.pos = translated_pos
        self.abs_pos = abs(self.pos)
        self.radius = radius
        self.semi_angle_size = Q_(atan(self.radius / self.abs_pos), "rad")

def compute_sphere_visibilities(spheres: list[Sphere]) -> dict[Sphere, float]:
    """Visibility from origin. 
    WARNING: do not consider multiple overlays!
    Visibility value will only be from the intersection of the largest two-spheres intersection"""
    visi = {s: 1 for s in spheres}
    for sphere in spheres:
        front_spheres = [s for s in spheres if s.abs_pos < sphere.abs_pos]
        for s in front_spheres:
            ang = s.pos.angle(sphere.pos)
            semi0 = sphere.semi_angle_size
            semi1 = s.semi_angle_size
            if ang >= semi0 + semi1:
                continue
            shaded_area = two_circles_intersection_area(semi0, semi1, ang)
            left_prop = 1 - shaded_area / (pi * semi0**2)
            visi[sphere] = left_prop if left_prop < visi[sphere] else visi[sphere]

