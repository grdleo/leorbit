from abc import ABC, abstractmethod
from typing import cast, overload

import numpy as np
import numpy.typing as npt

from leorbit.algorithms import sgp4
from leorbit.coordinates import Coordinates, OrbitalElements, Trajectory
from leorbit.frames import AbsoluteFrame
from leorbit.mathematics import Angle, Length, Quantity, Scalar, ScalarArray, Time, Vector3, Vector3Array, Velocity, normalize_angle
from leorbit.time import TimeInterval, Timestamp
from leorbit.utils import elements2orthogonal_gcrf, mean2true_anomaly

class Propagator(ABC):
    """Algorithm to propagate given orbital elements at given time"""
    
    def __init__(self, elements: OrbitalElements):
        self.elements = elements
    
    @overload
    def propagate(self, epoch: Timestamp) -> Coordinates:
        ...

    @overload
    def propagate(self, epoch: TimeInterval) -> Trajectory:
        ...
    
    @abstractmethod
    def propagate(self, epoch: TimeInterval | Timestamp) -> Trajectory | Coordinates:
        ...

class NoPropagator(Propagator):
    """A propagator that does not propagate, but always returns the same coordinates as given by the orbital elements"""

    @overload
    def propagate(self, epoch: Timestamp) -> Coordinates:
        ...

    @overload
    def propagate(self, epoch: TimeInterval) -> Trajectory:
        ...
    
    def propagate(self, epoch: TimeInterval | Timestamp) -> Trajectory | Coordinates:
        els = self.elements

        if isinstance(epoch, Timestamp):
            shift = (els.mean_motion * epoch.delta(els.epoch)).cast(Angle)
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
            time_line = epoch.to_time_stamps() - epoch.start.unixepoch * Quantity.second
            shift = (time_line * els.mean_motion).cast(Angle)
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
    def propagate(self, epoch: Timestamp) -> Coordinates:
        ...

    @overload
    def propagate(self, epoch: TimeInterval) -> Trajectory:
        ...
    
    def propagate(self, epoch: TimeInterval | Timestamp) -> Trajectory | Coordinates:
        tsince: npt.NDArray[np.float64]

        if isinstance(epoch, Timestamp):
            tsince = epoch.delta(self.elements.epoch).get_raw_array("minute")
        elif isinstance(epoch, TimeInterval):
            tsince = (epoch.to_time_stamps() - self.elements.epoch.unixepoch * Quantity.second).get_raw_array("minute")
        else:
            raise TypeError()

        output = sgp4(
            self.elements.compute_tuple,
            tsince.flatten()
        )

        if isinstance(epoch, Timestamp):
            x = float(np.asarray(output.x).reshape(-1)[0])
            y = float(np.asarray(output.y).reshape(-1)[0])
            z = float(np.asarray(output.z).reshape(-1)[0])
            vx = float(np.asarray(output.vx).reshape(-1)[0])
            vy = float(np.asarray(output.vy).reshape(-1)[0])
            vz = float(np.asarray(output.vz).reshape(-1)[0])

            return Coordinates(
                epoch,
                AbsoluteFrame.GCRF,
                cast(Vector3[Length], Vector3[Length].from_components(
                    x=x,
                    y=y,
                    z=z
                ).cast(Length)),
                cast(Vector3[Velocity], Vector3[Velocity].from_components(
                    x=vx,
                    y=vy,
                    z=vz
                ).cast(Velocity)),
            )
        elif isinstance(epoch, TimeInterval):
            return Trajectory(
                epoch,
                AbsoluteFrame.GCRF,
                cast(Vector3Array[Length], Vector3Array[Length].from_components(
                    x=output.x,
                    y=output.y,
                    z=output.z
                ).cast(Length)),
                cast(Vector3Array[Velocity], Vector3Array[Velocity].from_components(
                    x=output.vx,
                    y=output.vy,
                    z=output.vz
                ).cast(Velocity)),
            )