from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from functools import cached_property
from typing import Annotated, Any, Generic, Never, NewType, Self, TypeAlias, TypeVar, cast, overload
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

SomeDim = TypeVar("SomeDim", bound=Annotated[Dimension, DimensionObj])

TensorData = np.typing.NDArray[np.floating[Any]]

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
    
    def __eq__(self, o: object) -> bool:
        if not isinstance(o, DimensionalTensor):
            raise NotImplementedError()
        
        return (
            self._dimension == o._dimension 
            and self._values == o._values
        )
    
    def __neq__(self, o: DimensionalTensor) -> bool:
        return not self.__eq__(o)
    
    def __add__(self, o: DimensionalTensor | Number) -> DimensionalTensor:
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
    
    def __sub__(self, o: DimensionalTensor | Number) -> DimensionalTensor:
        o = ensure_tensor(o)
        
        if self._dimension != o._dimension:
            raise RuntimeError("Implement message")
        
        try:
            return self.__class__(
                self._values - o._values
            )
        except ValueError:
            raise RuntimeError("Implement message")
        
    def __rsub__(self, o: DimensionalTensor | Number) -> DimensionalTensor:
        o = ensure_tensor(o)
        
        if self._dimension != o._dimension:
            raise RuntimeError("Implement message")
        
        try:
            return self.__class__(
                o._values - self._values
            )
        except ValueError:
            raise RuntimeError("Implement message")
    
    def __mul__(self, o: DimensionalTensor | Number) -> DimensionalTensor:
        o = ensure_tensor(o)
        tensor_class = dimensional_tensor_class_factory(self.dimension * o.dimension)

        try:
            return tensor_class(
                self._values * o._values
            )
        except ValueError:
            raise RuntimeError("Implement message")
        
    def __rmul__(self, o: DimensionalTensor | Number) -> DimensionalTensor: # symetric
        return self.__mul__(o)
    
    def __truediv__(self, o: DimensionalTensor | Number) -> DimensionalTensor:
        o = ensure_tensor(o)
        tensor_class = dimensional_tensor_class_factory(self.dimension * o.dimension)

        try:
            return tensor_class(
                self._values / o._values
            )
        except ValueError:
            raise RuntimeError("Implement message")
        
    def __rtruediv__(self, o: DimensionalTensor | Number) -> DimensionalTensor:
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
            raise TypeError("Scalar must be instantiated with a dimension type: Scalar[LengthDim](value)")
        
        super().__init__(value)
    
    def __class_getitem__(cls, dim: SomeDim) -> type[Self]:
        return dim_to_tensor_class(dim, cls)

class Vector3(Generic[SomeDim], DimensionalTensor):
    def __init__(self, x: Number, y: Number, z: Number):
        if self._dimension is None:
            raise TypeError("Scalar must be instantiated with a dimension type: Vector3[LengthDim](value)")
        
        super().__init__(np.array([x, y, z]))
    
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
            raise TypeError("Scalar must be instantiated with a dimension type: Vector3[LengthDim](value)")
        
        super().__init__(np.array([[a, b, c], [d, e, f], [g, h, i]]))

    @property
    def inverse(self) -> Matrix33:
        raise NotImplementedError()
    
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

### TRANSFORMS ###

class Transform(Generic[DT1, DT2], ABC):
    @abstractmethod
    def do(self, tensor: DT1) -> DT2: ...

    @abstractmethod
    def undo(self, tensor: DT2) -> DT1: ...

    @abstractmethod
    def copy(self) -> Self: ...

class TransformIdentify(Generic[DT], Transform[DT, DT]):
    def do(self, tensor: DT) -> DT:
        return tensor
    
    def undo(self, tensor: DT) -> DT:
        return tensor
    
    def copy(self) -> Self:
        t = TransformIdentify[DT]()
        return cast(Self, t)

class TransformVector3Linear(Generic[SomeDim], Transform[Vector3[SomeDim], Vector3[SomeDim]]):
    def __init__(self, matrix: Matrix33[Dimensionless]):
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
    def __init__(self, matrix: Matrix33[Dimensionless], translation: Vector3[SomeDim]):
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
    def __init__(self, angle_rad: Number):
        c, s = np.cos(angle_rad), np.sin(angle_rad)
        rot_mat = Matrix33[Dimensionless](
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