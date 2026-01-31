from typing import Any
from leorbit2.mathematics import DimensionObj
import numpy as np

Number = int | float | np.floating[Any]
TensorData = np.typing.NDArray[np.floating[Any]]

class DimensionalTensor:
    _dimension: DimensionObj
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
        tensor_class = dimensional_tensor_class_factory(self._dimension * o._dimension)

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
        tensor_class = dimensional_tensor_class_factory(self._dimension * o._dimension)

        try:
            return tensor_class(
                self._values / o._values
            )
        except ValueError:
            raise RuntimeError("Implement message")
        
    def __rtruediv__(self, o: DimensionalTensor | Number) -> DimensionalTensor:
        o = ensure_tensor(o)
        tensor_class = dimensional_tensor_class_factory(self._dimension * o._dimension)

        try:
            return tensor_class(
                o._values / self._values
            )
        except ValueError:
            raise RuntimeError("Implement message")
        
    def __matmul__(self, o: DimensionalTensor) -> DimensionalTensor:
        o = ensure_tensor(o)
        tensor_class = dimensional_tensor_class_factory(self._dimension * o._dimension)

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

def dimensional_tensor_class_factory(dim: DimensionObj) -> type[DimensionalTensor]:
    # FIXME
    return type(
        "DimensionalTensor",
        (DimensionalTensor, ),
        dict(_dimension=dim)
    )