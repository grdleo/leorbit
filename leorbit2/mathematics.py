from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from functools import cached_property
from typing import Annotated, Any, Generic, Never, NewType, Self, TypeAlias, TypeIs, TypeVar, Union, cast, overload
import numpy as np
import numpy.typing

class Dimension:
    pass

Number: TypeAlias = int | float | np.floating[Any]

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
    
class Dim:
    DIM_TO_FACTORS: dict[DimensionObj, type[UnitRegistry]] = {}

    class UnitRegistry:
        pass

    class Dimensionless(UnitRegistry):
        _dim = DimensionObj()
    dimensionless = Annotated[Dimension, Dimensionless._dim]
    DIM_TO_FACTORS[Dimensionless._dim] = Dimensionless
    
    class Angle(UnitRegistry):
        _dim = DimensionObj()

        rad: float = 1 # base
        deg: float = np.pi / 180
    angle = Annotated[Dimension, Angle._dim]
    DIM_TO_FACTORS[Dimensionless._dim] = Angle
    
    class Length(UnitRegistry):
        _dim = DimensionObj(length=1)

        mm: float = 1e-3
        m: float = 1 # base
        km: float = 1e3
    length = Annotated[Dimension, Length._dim]
    DIM_TO_FACTORS[Length._dim] = Length
    
    class Time(UnitRegistry):
        _dim = DimensionObj(time=1)

        ms: float = 1e-3
        s: float = 1 # base
        min: float = 60 * s
        hour: float = 60 * min
        day: float = 24 * hour
        month: float = 30 * day
        year: float = 365 * day
    time = Annotated[Dimension, Time._dim]
    DIM_TO_FACTORS[Time._dim] = Time
    
    class Mass(UnitRegistry):
        _dim = DimensionObj(mass=1)

        g = 1e-3
        kg = 1 # base
        ton = 1e3
    mass = Annotated[Dimension, Mass._dim]
    DIM_TO_FACTORS[Mass._dim] = Mass

    class Velocity(UnitRegistry):
        _dim = DimensionObj(length=1, time=-1)

        m_s = 1 # m/s base 
        km_h = 1/3.6
    velocity = Annotated[Dimension, Velocity._dim]
    DIM_TO_FACTORS[Velocity._dim] = Velocity

    @staticmethod
    def get_factors(dobj: DimensionObj) -> type[UnitRegistry]:
        return Dim.DIM_TO_FACTORS[dobj]

SomeDim = TypeVar("SomeDim", bound=Annotated[Dimension, DimensionObj])
DimensionlessT = TypeVar("DimensionlessT", bound=Dim.dimensionless)
DimensionfullT = TypeVar("DimensionfullT", bound=Union[
    Dim.length, Dim.time, Dim.mass, Dim.velocity
])

TensorOrNumber = Union["DimensionalTensor", Number]

TensorData = np.typing.NDArray[np.floating[Any]]

def is_tensor_dim(tensor: DT, dim: type[Annotated[Dimension, DimensionObj]]) -> TypeIs[DT[SomeDim]]:
    dim_obj = cast(DimensionObj, dim.__metadata__[0])
    return tensor._dimension == dim_obj

class DimensionalTensor:
    _dimension: DimensionObj | None = None
    _values: TensorData

    def __init__(self, values: Number | TensorData):
        if isinstance(values, Number):
            self._values = np.array(values)
            return
        
        self._values = values

    @property
    def dimension(self) -> DimensionObj:
        """Access the dimension dynamically inside the class"""
        if self._dimension is None:
            raise RuntimeError("Dimension not set")
        return self._dimension
    
    @property
    def tensor_shape(self) -> tuple | tuple[int, ] | tuple[int, int]:
        return self._values.shape
    
    def copy(self) -> Self:
        return self.__class__(self._values.copy())
    
    def ensure_equal_dimensions(self, o: DimensionalTensor):
        if self._dimension != o._dimension:
            raise RuntimeError("Tensors' dimensions are not equal")
    
    def __eq__(self, o: object) -> bool:
        if not isinstance(o, DimensionalTensor):
            raise NotImplementedError()
        
        return (
            self._dimension == o._dimension 
            and self._values == o._values
        )
    
    def __neq__(self, o: DimensionalTensor) -> bool:
        return not self.__eq__(o)
    
    def __add__(self, o: TensorOrNumber) -> DimensionalTensor:
        o = ensure_tensor(o)
        
        if self._dimension != o._dimension:
            raise RuntimeError("Implement message")
        
        try:
            return self.__class__(
                self._values + o._values
            )
        except ValueError:
            raise RuntimeError("Implement message")
        
    def __radd__(self, o: Number) -> DimensionalTensor: # symetric
        return self.__add__(o)
    
    def __sub__(self, o: TensorOrNumber) -> DimensionalTensor:
        o = ensure_tensor(o)
        
        if self._dimension != o._dimension:
            raise RuntimeError("Implement message")
        
        try:
            return self.__class__(
                self._values - o._values
            )
        except ValueError:
            raise RuntimeError("Implement message")
        
    def __rsub__(self, o: TensorOrNumber) -> DimensionalTensor:
        o = ensure_tensor(o)
        
        if self._dimension != o._dimension:
            raise RuntimeError("Implement message")
        
        try:
            return self.__class__(
                o._values - self._values
            )
        except ValueError:
            raise RuntimeError("Implement message")
    
    def __mul__(self, o: TensorOrNumber) -> DimensionalTensor:
        o = ensure_tensor(o)
        tensor_class = dimensional_tensor_class_factory(self.dimension * o.dimension)

        try:
            return tensor_class(
                self._values * o._values
            )
        except ValueError:
            raise RuntimeError("Implement message")
        
    def __rmul__(self, o: TensorOrNumber) -> DimensionalTensor: # symetric
        return self.__mul__(o)
    
    def __truediv__(self, o: TensorOrNumber) -> DimensionalTensor:
        o = ensure_tensor(o)
        tensor_class = dimensional_tensor_class_factory(self.dimension * o.dimension)

        try:
            return tensor_class(
                self._values / o._values
            )
        except ValueError:
            raise RuntimeError("Implement message")
        
    def __rtruediv__(self, o: TensorOrNumber) -> DimensionalTensor:
        o = ensure_tensor(o)
        tensor_class = dimensional_tensor_class_factory(self.dimension * o.dimension)

        try:
            return tensor_class(
                o._values / self._values
            )
        except ValueError:
            raise RuntimeError("Implement message")
        
    def __matmul__(self, o: DimensionalTensor) -> DimensionalTensor:
        o = ensure_tensor(o)
        tensor_class = dimensional_tensor_class_factory(self.dimension * o.dimension)

        try:
            return tensor_class(
                o._values @ self._values
            )
        except ValueError:
            raise RuntimeError("Implement message")

def ensure_tensor(el: Number | TensorData | DimensionalTensor) -> DimensionalTensor:
    """ensure tensor. if not a tensor object, creates a dimless tensor"""

    if isinstance(el, DimensionalTensor):
        return el
    elif isinstance(el, Number | TensorData):
        dimensionless = DimensionObj()
        return dimensional_tensor_class_factory(dimensionless)(el)
    
    raise RuntimeError("...")

DT = TypeVar("DT", bound=DimensionalTensor)
DT1 = TypeVar("DT1", bound=DimensionalTensor)
DT2 = TypeVar("DT2", bound=DimensionalTensor)

def dimensional_tensor_class_factory(
    dim: DimensionObj, 
    parent_class: type[DimensionalTensor] = DimensionalTensor
) -> type[DimensionalTensor]:
    # FIXME
    return type(
        "DimensionalTensor",
        (parent_class, ),
        dict(_dimension=dim)
    )

def dim_to_tensor_class(dim: SomeDim, tensor_class: type[DT]) -> type[DT]:
     # Extract the DimensionObj from the Annotated type
    if hasattr(dim, '__metadata__'):
        dimension_obj = dim.__metadata__[0]
    else:
        raise TypeError(f"Dimension type must be Annotated with a DimensionObj")
    
    # FIXME
    return type(
        "DimensionalTensor",
        (tensor_class, ),
        dict(_dimension=dim)
    )

class Scalar(Generic[SomeDim], DimensionalTensor):
    def __init__(self, value: Number):
        if self._dimension is None:
            raise TypeError("Scalar must be instantiated with a dimension type: Scalar[Dim.length](value)")
        
        super().__init__(value)

    @property
    def base_units_value(self) -> Number:
        return np.float64(self._values)
    
    def magnitude(self, unit: str = "1") -> Number:
        if self._dimension is None:
            raise TypeError("Scalar must be instantiated with a dimension type: Scalar[Dim.length](value)")
        
        base_value = np.float64(self._values)

        if unit == "1":
            return base_value
        
        dim_factors = Dim.get_factors(self._dimension)
        factor = getattr(dim_factors, unit, None)
        if factor is None:
            raise ValueError(f"No unit '{unit}' for dimension '{self._dimension}'")
        
        return base_value / factor

    @overload
    def __add__(self, o: Number) -> Scalar[DimensionlessT] | Never:
        if isinstance(self, Scalar[Dim.dimensionless]):
            return Scalar[Dim.dimensionless](0)
        return None

    @overload
    def __add__(self: Scalar[DimensionfullT], o: Number) -> Never: ...

    @overload
    def __add__(self, o: DT) -> DT: ...

    def __add__(self, o: TensorOrNumber) -> DimensionalTensor:
        return super().__add__(o)
    
    @overload
    def __sub__(self, o: DT) -> DT: ...
    
    @overload
    def __sub__(self, o: Number) -> Scalar[SomeDim]: ...

    def __sub__(self, o: TensorOrNumber) -> DimensionalTensor:
        return super().__sub__(o)
    
    @overload
    def __mul__(self, o: DT) -> DT: ...
    
    @overload
    def __mul__(self, o: Number) -> Scalar[SomeDim]: ...

    def __mul__(self, o: TensorOrNumber) -> DimensionalTensor:
        return super().__mul__(o)
    
    @overload
    def __truediv__(self, o: DT) -> Never: ...
    
    @overload
    def __truediv__(self, o: Number | Scalar[SomeDim]) -> Scalar[SomeDim]: ...

    def __truediv__(self, o: TensorOrNumber) -> DimensionalTensor:
        return super().__truediv__(o)
    
    def __lt__(self, o: Scalar[SomeDim]) -> bool:
        self.ensure_equal_dimensions(o)
        return bool(self._values < o._values)

    def __le__(self, o: Scalar[SomeDim]) -> bool:
        self.ensure_equal_dimensions(o)
        return bool(self._values <= o._values)

    def __gt__(self, o: Scalar[SomeDim]) -> bool:
        self.ensure_equal_dimensions(o)
        return bool(self._values > o._values)

    def __ge__(self, o: Scalar[SomeDim]) -> bool:
        self.ensure_equal_dimensions(o)
        return bool(self._values >= o._values)
    
    def __class_getitem__(cls, dim: SomeDim) -> type[Self]:
        return dim_to_tensor_class(dim, cls)

class Vector3(Generic[SomeDim], DimensionalTensor):
    def __init__(self, x: Number, y: Number, z: Number):
        if self._dimension is None:
            raise TypeError("Scalar must be instantiated with a dimension type: Vector3[Dim.length](value)")
        
        super().__init__(np.array([x, y, z]).reshape((3,1)))
    
    def __class_getitem__(cls, dim: SomeDim) -> type[Self]:
        return dim_to_tensor_class(dim, cls)
    
class Matrix33(Generic[SomeDim], DimensionalTensor):
    def __init__(self,
        a: Number, b: Number, c: Number,
        d: Number, e: Number, f: Number,
        g: Number, h: Number, i: Number,
    ):
        """order: by lines"""
        if self._dimension is None:
            raise TypeError("Scalar must be instantiated with a dimension type: Vector3[Dim.length](value)")
        
        super().__init__(np.array([[a, b, c], [d, e, f], [g, h, i]]))

    @property
    def inverse(self) -> Matrix33:
        raise NotImplementedError()
    
############################################################
    
class QuantityMeta(type):
    def __getattr__(cls, name: str) -> Scalar:
        if name == "rad":
            return Scalar[Dim.dimensionless](1)
        elif name == "deg":
            return Scalar[Dim.dimensionless](Dim.Angle.deg)
        elif name == "m":
            return Scalar[Dim.length](1)
        elif name == "km":
            return Scalar[Dim.length](Dim.Length.km)
        elif name == "s":
            return Scalar[Dim.time](1)
        elif name == "min":
            return Scalar[Dim.time](Dim.Time.min)
        elif name == "hour":
            return Scalar[Dim.time](Dim.Time.hour)
        elif name == "day":
            return Scalar[Dim.time](Dim.Time.day)
        elif name == "month":
            return Scalar[Dim.time](Dim.Time.month)
        elif name == "year":
            return Scalar[Dim.time](Dim.Time.year)
        
        raise ValueError(f"No unit named '{name}'")
        

class Quantity(metaclass=QuantityMeta):
    # ANGLES

    rad: Scalar[Dim.dimensionless]
    """radians"""

    deg: Scalar[Dim.dimensionless]
    """degrees"""

    # DISTANCES

    m: Scalar[Dim.length]
    """meter"""

    km: Scalar[Dim.length]
    """kilometer"""

    # DURATIONS

    s: Scalar[Dim.time]
    """second"""

    min: Scalar[Dim.time]
    """minute"""

    hour: Scalar[Dim.time]
    """hour"""

    day: Scalar[Dim.time]
    """day"""

    month: Scalar[Dim.time]
    """month (30 days)"""

    year: Scalar[Dim.time]
    """year (365 days)"""

def quantity(nb_and_unit: str) -> Scalar:
    nb, unit_name = nb_and_unit.split(" ")
    unit: Scalar = getattr(Quantity, unit_name)
    return cast(Scalar, unit * float(nb))


### TRANSFORMS ###

class Transform(Generic[DT1, DT2], ABC):
    @abstractmethod
    def do(self, tensor: DT1) -> DT2: ...

    @abstractmethod
    def undo(self, tensor: DT2) -> DT1: ...

    @abstractmethod
    def copy(self) -> Self: ...

    def reverse(self) -> Transform[DT2, DT1]:
        t = self.copy()

        do = t.do
        undo = t.undo

        tt = cast(Transform[DT2, DT1], t)

        tt.do = undo  # type: ignore
        tt.undo = do  # type: ignore

        return tt

class TransformIdentify(Generic[DT], Transform[DT, DT]):
    def do(self, tensor: DT) -> DT:
        return tensor
    
    def undo(self, tensor: DT) -> DT:
        return tensor
    
    def copy(self) -> Self:
        t = TransformIdentify[DT]()
        return cast(Self, t)

class TransformVector3Linear(Generic[SomeDim], Transform[Vector3[SomeDim], Vector3[SomeDim]]):
    def __init__(self, matrix: Matrix33[Dim.dimensionless]):
        self.matrix = matrix
    
    def do(self, v: Vector3[SomeDim]) -> Vector3[SomeDim]:
        return cast(Vector3[SomeDim], self.matrix * v)
    
    def undo(self, v: Vector3[SomeDim]) -> Vector3[SomeDim]:
        return cast(Vector3[SomeDim], self.matrix.inverse * v)
    
    def copy(self) -> Self:
        t = TransformVector3Linear[SomeDim](
            matrix=self.matrix.copy()
        )
        return cast(Self, t)

class TransformVector3Affine(Generic[SomeDim], Transform[Vector3[SomeDim], Vector3[SomeDim]]):
    def __init__(self, matrix: Matrix33[Dim.dimensionless], translation: Vector3[SomeDim]):
        self.matrix = matrix
        self.translation = translation
    
    def do(self, v: Vector3[SomeDim]) -> Vector3[SomeDim]:
        result = self.matrix * v + self.translation
        return cast(Vector3[SomeDim], result)
    
    def undo(self, v: Vector3[SomeDim]) -> Vector3[SomeDim]:
        result = self.matrix.inverse * (v - self.translation)
        return cast(Vector3[SomeDim], result)
    
    def copy(self) -> Self:
        t = TransformVector3Affine[SomeDim](
            matrix=self.matrix.copy(), 
            translation=self.translation.copy()
        )
        return cast(Self, t)
    
class TransformVector3RotationZ(Generic[SomeDim], TransformVector3Linear[SomeDim]):
    def __init__(self, angle_rad: Number | Scalar[Dim.angle]):
        if isinstance(angle_rad, Scalar):
            angle_rad = angle_rad.magnitude("rad")
        
        c, s = np.cos(angle_rad), np.sin(angle_rad)
        rot_mat = Matrix33[Dim.dimensionless](
            c, -s, 0,
            s, c, 0,
            0, 0, 1
        )

        super().__init__(rot_mat) 
    
class TransformChain(Generic[DT1, DT2], Transform[DT1, DT2]):
    def __init__(self, *transforms: Transform[DimensionalTensor, DimensionalTensor]):
        """executed in the given order"""
        self.transforms = list(transforms)

    def do(self, tensor: DT1) -> DT2:
        result = cast(DimensionalTensor, tensor)
        for t in self.transforms:
            result = t.do(result)
        return cast(DT2, result)
    
    def undo(self, tensor: DT2) -> DT1:
        result = cast(DimensionalTensor, tensor)
        for t in reversed(self.transforms):
            result = t.undo(result)
        return cast(DT1, result)
    
    def copy(self) -> Self:
        t = TransformChain[DT1, DT2](
            *(t.copy() for t in self.transforms)
        )

        return cast(Self, t)
    




a = Scalar[Dim.mass](123)
b = Scalar[Dim.mass](321)
c = b + 1