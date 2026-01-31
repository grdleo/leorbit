from dataclasses import dataclass, field
from typing import Annotated, Any, Generic, NewType, Self, TypeAlias, TypeVar, overload
import numpy as np
import numpy.typing

class Dimension:
    pass

Number: TypeAlias = int | float

@dataclass(eq=True, frozen=True)
class DimensionObj:
    length: int = 0
    time: int = 0
    mass: int = 0

    def __repr__(self) -> str:
        components = [
            ("" if self.length == 0 else ("m" if self.length == 1 else f"m^{self.length}")),
            ("" if self.time == 0 else ("s" if self.time == 1 else f"s^{self.time}")),
            ("" if self.mass == 0 else ("kg" if self.mass == 1 else f"kg^{self.mass}"))
        ]

        return " ".join(components)

    @property
    def dimensionless(self) -> bool:
        return (
            self.length
            == self.time
            == self.mass
            == 0
        )
    
    def __mul__(self, o: DimensionObj) -> DimensionObj:
        return DimensionObj(
            length=self.length + o.length,
            time=self.time + o.time,
            mass=self.mass + o.mass
        )
    
    def __truediv__(self, o: DimensionObj) -> DimensionObj:
        return DimensionObj(
            length=self.length - o.length,
            time=self.time - o.time,
            mass=self.mass - o.mass
        )
    
@dataclass
class DimensionRegister:
    dimensionless: DimensionObj = field(default=DimensionObj(), init=False)
    length: DimensionObj = field(default=DimensionObj(length=1), init=False)
    time: DimensionObj = field(default=DimensionObj(time=1), init=False)
    mass: DimensionObj = field(default=DimensionObj(mass=1), init=False)
    velocity: DimensionObj = field(default=DimensionObj(length=1, time=-1), init=False)
    
REGISTER = DimensionRegister()

Dimensionless = Annotated[Dimension, REGISTER.dimensionless]
LengthDim = Annotated[Dimension, REGISTER.length]
TimeDim = Annotated[Dimension, REGISTER.time]
MassDim = Annotated[Dimension, REGISTER.mass]
VelocityDim = Annotated[Dimension, REGISTER.velocity]

AnyDim = TypeVar("AnyDim", bound=Annotated[Dimension, DimensionObj])

class S:
    _dimension: DimensionObj
    _value: Number

    def __init__(self, value: Number):
        self._value = value

    @property
    def dimension(self) -> DimensionObj:
        """Access the dimension dynamically inside the class"""
        if self._dimension is None:
            raise RuntimeError("Dimension not set")
        return self._dimension
    
    def __repr__(self) -> str:
        return f"{self._value} {self._dimension}"
    
    def __eq__(self, o: object) -> bool:
        if not isinstance(o, S):
            raise NotImplementedError()
        
        return (
            self._dimension == o._dimension 
            and self._value == o._value
        )
    
    def __neq__(self, o: Self) -> bool:
        return not self.__eq__(o)
    
    def __mul__(self, o: S | Number) -> S:
        if isinstance(o, Number):
            return self.__class__(self._value * o)
        elif isinstance(o, S):
            return scalar_class_factory(
                self._dimension * o._dimension
            )(self._value * o._value)
        
        raise TypeError()
    
    def __rmul__(self, o: S | Number) -> S:
        return self.__mul__(o)
    
    def __add__(self, o: S | Number) -> S:
        if isinstance(o, Number):
            if not self._dimension.dimensionless:
                raise TypeError()
            return self.__class__(self._value + o)
        if self._dimension != o._dimension:
            raise TypeError()
        
        return self.__class__(self._value + o._value)
    
    def __radd__(self, o: S | Number) -> S:
        return self.__add__(o)
    
    def __sub__(self, o: S | Number) -> S:
        if isinstance(o, Number):
            if not self._dimension.dimensionless:
                raise TypeError()
            return self.__class__(self._value - o)
        if self._dimension != o._dimension:
            raise TypeError()
        
        return self.__class__(self._value - o._value)
    
    def __rsub__(self, o: S | Number) -> S:
        if isinstance(o, Number):
            if not self._dimension.dimensionless:
                raise TypeError()
            return self.__class__(o - self._value)
        if self._dimension != o._dimension:
            raise TypeError()
        
        return self.__class__(o._value - self._value)
    
    def __truediv__(self, o: S | Number) -> S:
        if isinstance(o, Number):
            return self.__class__(self._value / o)
        elif isinstance(o, S):
            return scalar_class_factory(
                self._dimension / o._dimension
            )(self._value / o._value)
        
        raise TypeError()
    
    def __rtruediv__(self, numerator: Number) -> S:
        return scalar_class_factory(
                REGISTER.dimensionless / self._dimension
            )(numerator / self._value)
        
def scalar_class_factory(dim: DimensionObj) -> type[S]:
    # FIXME cache the types
    return type(
        "S", 
        (S, ),
        dict(_dimension=dim)
    )

class Scalar(Generic[AnyDim], S):
    def __init__(self, *args, **kwargs):
        raise TypeError("Scalar must be instantiated with a dimension type: Scalar[LengthDim](value)")
    
    def __class_getitem__(cls, dim: type) -> type:
        """Capture the dimension type when user does Scalar[SomeDim]"""
        # Extract the DimensionObj from the Annotated type
        if hasattr(dim, '__metadata__'):
            dimension_obj = dim.__metadata__[0]
        else:
            raise TypeError(f"Dimension type must be Annotated with a DimensionObj")
        
        return scalar_class_factory(dimension_obj)
    
class QuantityMeta(type):
    def __getattr__(cls, name: str) -> Any:
        if name == "m":
            return Scalar[LengthDim](1)
        elif name == "km":
            return Scalar[LengthDim](1000)

class Quantity(metaclass=QuantityMeta):
    m: Scalar[LengthDim]
    """meter"""

    km: Scalar[LengthDim]
    """kilometer"""

class AS:
    _dimension: DimensionObj
    _values: np.typing.NDArray[np.floating[Any]]

    def __init__(self, values: np.typing.NDArray[np.floating[Any]]):
        self._values = values

    @property
    def dimension(self) -> DimensionObj:
        """Access the dimension dynamically inside the class"""
        if self._dimension is None:
            raise RuntimeError("Dimension not set")
        return self._dimension
    
    def __repr__(self) -> str:
        return f"{self._values} {self._dimension}"
    
    def __eq__(self, o: object) -> bool:
        if not isinstance(o, AS):
            raise NotImplementedError()
        
        return (
            self._dimension == o._dimension 
            and self._values == o._values
        )
    
    def __neq__(self, o: Self) -> bool:
        return not self.__eq__(o)
    
    def __mul__(self, o: AS | S | Number) -> AS:
        if isinstance(o, Number):
            return self.__class__(self._values * o)
        elif isinstance(o, S):
            return array_scalar_class_factory(
                self._dimension * o._dimension
            )(self._values * o._value)
        elif isinstance(o, AS):
            return array_scalar_class_factory(
                self._dimension * o._dimension
            )(self._values * o._values)
        
        raise TypeError()
    
    def __rmul__(self, o: AS | Number) -> AS:
        return self.__mul__(o)
    
    def __add__(self, o: AS | S | Number) -> AS:
        if isinstance(o, Number):
            if not self._dimension.dimensionless:
                raise TypeError()
            return self.__class__(self._values + o)
        elif isinstance(o, S):
            return self.__class__(self._values + o._value)
        if self._dimension != o._dimension:
            raise TypeError()
        
        return self.__class__(self._values + o._values)
    
    def __radd__(self, o: AS | Number) -> AS:
        return self.__add__(o)
    
    def __sub__(self, o: AS | S | Number) -> AS:
        if isinstance(o, Number):
            if not self._dimension.dimensionless:
                raise TypeError()
            return self.__class__(self._values - o)
        elif isinstance(o, S):
            return self.__class__(self._values - o._value)
        if self._dimension != o._dimension:
            raise TypeError()
        
        return self.__class__(self._values - o._values)
    
    def __rsub__(self, o: AS | S | Number) -> AS:
        if isinstance(o, Number):
            if not self._dimension.dimensionless:
                raise TypeError()
            return self.__class__(o - self._values)
        elif isinstance(o, S):
            return self.__class__(o._value - self._values)
        if self._dimension != o._dimension:
            raise TypeError()
        
        return self.__class__(o._values - self._values)
    
    def __truediv__(self, o: AS | S | Number) -> AS:
        if isinstance(o, Number):
            return self.__class__(self._values / o)
        elif isinstance(o, S):
            return array_scalar_class_factory(
                self._dimension / o._dimension
            )(self._values / o._value)
        elif isinstance(o, AS):
            return array_scalar_class_factory(
                self._dimension / o._dimension
            )(self._values / o._values)
        
        raise TypeError()
    
    def __rtruediv__(self, numerator: Number) -> S:
        return array_scalar_class_factory(
                REGISTER.dimensionless / self._dimension
            )(numerator / self._values)
    
def get_tensor_value(tensor: Number | S | AS) -> Number | np.typing.NDArray[np.floating]:
    if isinstance(tensor, Number):
        return np.float64(tensor)
    elif isinstance(tensor, S):
        return tensor._value
    elif isinstance(tensor, AS):
        return tensor._values
    
    raise TypeError()

    
a = Scalar[TimeDim](5.0)
b = Scalar[MassDim](3.2)
c = 1 / a

print(c)