from abc import ABC, abstractmethod
from dataclasses import dataclass
from fractions import Fraction
from pyclbr import Class
from typing import Any, ClassVar, Generic, Literal, NamedTuple, Never, Self, TypeAlias, TypeGuard, TypeIs, TypeVar, cast, overload

import numpy as np

class DimEls:
    def __init__(self,
        length: Fraction | int = 0,
        time: Fraction | int = 0,
        mass: Fraction | int = 0,
    ):
        self.__length = Fraction(length)
        self.__time = Fraction(time)
        self.__mass = Fraction(mass)

    @property
    def length(self) -> Fraction:
        return self.__length
    
    @property
    def time(self) -> Fraction:
        return self.__time
    
    @property
    def mass(self) -> Fraction:
        return self.__mass

    @property
    def dimensionless(self) -> bool:
        return (
            0
            == self.length
            == self.time
            == self.mass
        )
    
    def __eq__(self, o: object) -> bool:
        if not isinstance(o, DimEls):
            return False
        
        return (
            self.length == o.length
            and self.time == o.time
            and self.mass == o.mass
        )
    
    def __mul__(self, o: DimEls) -> DimEls:
        return DimEls(
            length=self.length + o.length,
            time=self.time + o.time,
            mass=self.mass + o.mass
        )
    
    def __truediv__(self, o: DimEls) -> DimEls:
        return DimEls(
            length=self.length - o.length,
            time=self.time - o.time,
            mass=self.mass - o.mass
        )
    
    def __pow__(self, p: Fraction | int | float) -> DimEls:
        p = Fraction(p)

        return DimEls(
            length=self.length * p,
            time=self.time * p,
            mass=self.mass * p
        )

class Dim:
    """Base class for dimensions"""
    _d: ClassVar[DimEls | None] = None


class D:
    """Dimensions registry."""

    class Dimless(Dim):
        """1"""
        _d = DimEls()

    class Angle(Dim):
        """rad"""
        _d = DimEls()

        rad: ClassVar[float] = 1 # base
        deg: ClassVar[float] = np.pi / 180

    class AngularVelocity(Dim):
        """rad/s"""
        _d = DimEls(time=-1)

    class AngularAcc(Dim):
        """rad/s**2"""
        _d = DimEls(time=-2)

    class AngularJerk(Dim):
        """rad/s**2"""
        _d = DimEls(time=-3)

    class Length(Dim):
        """m"""
        _d = DimEls(length=1)

        meter: ClassVar[float] = 1 # base

        milli_meter: ClassVar[float] = 1e-3 * meter
        kilo_meter: ClassVar[float] = 1e3 * meter

    class InvLength(Dim):
        """m**-1"""
        _d = DimEls(length=-1)

    class Time(Dim):
        """s"""
        _d = DimEls(time=1)

        second: ClassVar[float] = 1 # base

        milli_second: ClassVar[float] = 1e-3 * second
        minute: ClassVar[float] = 60 * second
        hour: ClassVar[float] = 60 * minute
        day: ClassVar[float] = 24 * hour
        month: ClassVar[float] = 30 * day
        year: ClassVar[float] = 365 * day

    class Mass(Dim):
        """kg"""
        _d = DimEls(mass=1)

        kilo_gram: ClassVar[float] = 1 # base

        gram: ClassVar[float] = 1e-3 * kilo_gram
        milli_gram: ClassVar[float] = 1e-3 * gram
        metric_ton: ClassVar[float] = 1e3 * kilo_gram

    class Velocity(Dim):
        """m.s**-1"""
        _d = DimEls(length=1, time=-1)

        meter_per_second: ClassVar[float] = 1 # base

        kilo_meter_per_hour: ClassVar[float] = meter_per_second / 3.6

    class Acceleration(Dim):
        """m.s**-2"""
        _d = DimEls(length=1, time=-2)

    class Force(Dim):
        """kg.m.s**-2"""
        _d = DimEls(mass=1, length=1, time=-2)

        newton: ClassVar[float] = 1 # base

    class GrativationnalParam(Dim):
        """m**3.s**-2"""
        _d = DimEls(length=3, time=-2)

    @staticmethod
    def registered_dimensions() -> dict[DimEls, type[Dim]]:
        return {
            dim_cls._d: dim_cls
            for dim_cls in D.__dict__.values()
            if issubclass(dim_cls, Dim)
        }

Number = float | int | np.floating
SomeDim = TypeVar("SomeDim", bound=Dim)
SomeOtherDim = TypeVar("SomeOtherDim", bound=Dim)
SomeDimFull = TypeVar("SomeDimFull", bound=D.Length | D.Time | D.Velocity)
TensorData = np.typing.NDArray[np.floating[Any]]

class ProductDim[SomeDim, SomeOtherDim](Dim):
    ...

class QuotientDim[SomeDim, SomeOtherDim](Dim):
    ...
    
class PowerDim[SomeDim, Numerator, Denominator](Dim):
    ...


N3: TypeAlias = Literal[-3]
N2: TypeAlias = Literal[-2]
N1: TypeAlias = Literal[-1]
OO: TypeAlias = Literal[0]
P1: TypeAlias = Literal[1]
P2: TypeAlias = Literal[2]
P3: TypeAlias = Literal[3]