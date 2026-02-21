from abc import ABC, abstractmethod
from typing import cast, overload

import numpy as np

from leorbit2.algorithms import TimeInterval, sgp4
from leorbit2.coordinates import Coordinates, OrbitalElements, Trajectory
from leorbit2.frames import AbsoluteFrame
from leorbit2.m import D, Scalar, ScalarArray, Vector3, Vector3Array, normalize_angle
from leorbit2.time import Time
from leorbit2.utils import elements2orthogonal_gcrf, mean2true_anomaly

import numpy.typing as npt

class Propagator(ABC):
    """Algorithm to propagate given orbital elements at given time"""
    
    def __init__(self, elements: OrbitalElements):
        self.elements = elements
    
    @overload
    def propagate(self, epoch: Time) -> Coordinates:
        ...

    @overload
    def propagate(self, epoch: TimeInterval) -> Trajectory:
        ...
    
    @abstractmethod
    def propagate(self, epoch: TimeInterval | Time) -> Trajectory | Coordinates:
        ...

class NoPropagator(Propagator):
    """A propagator that does not propagate, but always returns the same coordinates as given by the orbital elements"""

    @overload
    def propagate(self, epoch: Time) -> Coordinates:
        ...

    @overload
    def propagate(self, epoch: TimeInterval) -> Trajectory:
        ...
    
    def propagate(self, epoch: TimeInterval | Time) -> Trajectory | Coordinates:
        els = self.elements

        if isinstance(epoch, Time):
            shift = (els.mean_motion * epoch.delta(els.epoch)).cast(D.Angle)
            shifted_M0 = normalize_angle(els.mean_anomaly + shift)
            shifted_nu = mean2true_anomaly(els.eccentricity, shifted_M0)
            pos, vel = elements2orthogonal_gcrf(
                shifted_nu,
                els.eccentricity,
                els.semi_major_axis,
                els.ra_of_asc_node,
                els.arg_of_pericenter,
                els.inclination
            )
            return Coordinates(epoch, AbsoluteFrame.GCRF, pos, vel)
        
        elif isinstance(epoch, TimeInterval):
            time_line = epoch.to_time_stamps() - Scalar[D.Time](epoch.start.unixepoch)
            shift = (time_line * els.mean_motion).cast(D.Angle)
            shifted_M0 = normalize_angle(els.mean_anomaly + shift)
            shifted_nu = mean2true_anomaly(els.eccentricity, shifted_M0)
            pos, vel = elements2orthogonal_gcrf(
                shifted_nu,
                els.eccentricity,
                els.semi_major_axis,
                els.ra_of_asc_node,
                els.arg_of_pericenter,
                els.inclination
            )
            return Trajectory(epoch, AbsoluteFrame.GCRF, pos, vel)
        
        else:
            raise TypeError()
        
class SGP4(Propagator):
    """A propagator that uses the SGP4 algorithm to propagate the orbital elements. 
    Note: SGP4 only works for Earth satellites, so the absolute frame of the returned coordinates is always GCRF"""

    @overload
    def propagate(self, epoch: Time) -> Coordinates:
        ...

    @overload
    def propagate(self, epoch: TimeInterval) -> Trajectory:
        ...
    
    def propagate(self, epoch: TimeInterval | Time) -> Trajectory | Coordinates:
        tsince: npt.NDArray[np.float64]

        if isinstance(epoch, Time):
            tsince = epoch.delta(self.elements.epoch).get_raw_array("minute")
        elif isinstance(epoch, TimeInterval):
            tsince = (epoch.to_time_stamps() - Scalar[D.Time](self.elements.epoch.unixepoch)).get_raw_array("minute")
        else:
            raise TypeError()

        output = sgp4(
            self.elements.compute_tuple,
            tsince.flatten()
        )

        if isinstance(epoch, Time):
            x = float(np.asarray(output.x).reshape(-1)[0])
            y = float(np.asarray(output.y).reshape(-1)[0])
            z = float(np.asarray(output.z).reshape(-1)[0])
            vx = float(np.asarray(output.vx).reshape(-1)[0])
            vy = float(np.asarray(output.vy).reshape(-1)[0])
            vz = float(np.asarray(output.vz).reshape(-1)[0])

            return Coordinates(
                epoch,
                AbsoluteFrame.GCRF,
                Vector3[D.Length].from_components(
                    x=x,
                    y=y,
                    z=z
                ).cast(D.Length),
                Vector3[D.Velocity].from_components(
                    x=vx,
                    y=vy,
                    z=vz
                ).cast(D.Velocity),
            )
        elif isinstance(epoch, TimeInterval):
            return Trajectory(
                epoch,
                AbsoluteFrame.GCRF,
                Vector3Array[D.Length].from_components(
                    x=output.x,
                    y=output.y,
                    z=output.z
                ).cast(D.Length),
                Vector3Array[D.Velocity].from_components(
                    x=output.vx,
                    y=output.vy,
                    z=output.vz
                ).cast(D.Velocity),
            )