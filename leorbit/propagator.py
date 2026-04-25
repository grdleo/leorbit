from abc import ABC, abstractmethod
from typing import overload

import numpy as np
import numpy.typing as npt

from leorbit.algorithms import sgp4
from leorbit.coordinates import Coordinates, OrbitalElements, Trajectory
from leorbit.frames import AbsoluteFrame
from leorbit.mathematics import Angle, Length, Quantity, Tensor, TensorBound, TensorKind, Time, normalize_angle
from leorbit.time import TimeInterval, Timestamp
from leorbit.utils import elements2orthogonal_gcrf, mean2true_anomaly

class Propagator(ABC):
    """Algorithm to propagate given orbital elements at given time"""
    
    def __init__(self, elements: OrbitalElements):
        """Store orbital elements used by this propagator implementation."""
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
        """Propagate by analytical Keplerian conversion without perturbations."""
        els = self.elements

        if isinstance(epoch, Timestamp):
            shift = (els.mean_motion * epoch.delta(els.epoch)).secure(Angle, TensorKind.SCALAR)
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
            shift = (time_line * els.mean_motion).secure(Angle, TensorKind.SCALAR)
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
        """Propagate state using the SGP4 numerical model."""
        tsince: npt.NDArray[np.float64]

        if isinstance(epoch, Timestamp):
            tsince = np.asarray([epoch.delta(self.elements.epoch).scalar.value("minute")], dtype=np.float64)
        elif isinstance(epoch, TimeInterval):
            tsince = (epoch.to_time_stamps() - self.elements.epoch.unixepoch * Quantity.second).raw_data_array("minute")
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
                Tensor(np.asarray([[x], [y], [z]], dtype=np.float64), Length),
                Tensor(np.asarray([[vx], [vy], [vz]], dtype=np.float64), Length / Time),
            )
        elif isinstance(epoch, TimeInterval):
            return Trajectory(
                epoch,
                AbsoluteFrame.GCRF,
                Tensor(np.asarray([output.x, output.y, output.z], dtype=np.float64), Length),
                Tensor(np.asarray([output.vx, output.vy, output.vz], dtype=np.float64), Length / Time),
            )