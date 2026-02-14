import math
from typing import Any, TypeVar, overload, Literal

import numpy as np
from leorbit2.mathematics.dimensions import D, Dim, Number, SomeDim, PowerDim, P1, P2
from leorbit2.mathematics.quantity import Quantity
from leorbit2.mathematics.scalar import Scalar
from leorbit2.mathematics.scalar_array import ScalarArray
from leorbit2.mathematics.tensor import Tensor

TWELF_PI = math.pi / 12
TWOPI = 2 * math.pi

FULL_REV = (2 * math.pi) * Quantity.rad
HALF_REV = FULL_REV / 2

@overload
def square(tensor: Scalar[D.Dimless]) -> Scalar[D.Dimless]: ... # type: ignore

@overload
def square(tensor: Scalar[PowerDim[SomeDim, P1, P2]]) -> Scalar[SomeDim]: ... # type: ignore

@overload
def square(tensor: Scalar[SomeDim]) -> Scalar[PowerDim[SomeDim, P2, P1]]: ... # type: ignore

@overload
def square(tensor: ScalarArray[D.Dimless]) -> ScalarArray[D.Dimless]: ... # type: ignore

@overload
def square(tensor: ScalarArray[PowerDim[SomeDim, P1, P2]]) -> ScalarArray[SomeDim]: ... # type: ignore

@overload
def square(tensor: ScalarArray[SomeDim]) -> ScalarArray[PowerDim[SomeDim, P2, P1]]: ... # type: ignore

def square(tensor: Tensor[Any]) -> Tensor[Any]:
    return tensor._base_tensor_class.dimensionalize(
        tensor.dim_coords ** 2
    )(
        np.square(tensor._values)
    )

@overload
def sqrt(tensor: Scalar[D.Dimless]) -> Scalar[D.Dimless]: ... # type: ignore

@overload
def sqrt(tensor: Scalar[PowerDim[SomeDim, P2, P1]]) -> Scalar[SomeDim]: ... # type: ignore

@overload
def sqrt(tensor: Scalar[SomeDim]) -> Scalar[PowerDim[SomeDim, P1, P2]]: ... # type: ignore

@overload
def sqrt(tensor: ScalarArray[D.Dimless]) -> ScalarArray[D.Dimless]: ... # type: ignore

@overload
def sqrt(tensor: ScalarArray[PowerDim[SomeDim, P2, P1]]) -> ScalarArray[SomeDim]: ... # type: ignore

@overload
def sqrt(tensor: ScalarArray[SomeDim]) -> ScalarArray[PowerDim[SomeDim, P1, P2]]: ... # type: ignore

def sqrt(tensor: Tensor[Any]) -> Tensor[Any]:
    return tensor._base_tensor_class.dimensionalize(
        tensor.dim_coords ** .5
    )(
        np.sqrt(tensor._values)
    )

@overload
def cos(tensor: Scalar[D.Angle]) -> Scalar[D.Dimless]: ... # type: ignore

@overload
def cos(tensor: ScalarArray[D.Angle]) -> ScalarArray[D.Dimless]: ... # type: ignore

def cos(tensor: Tensor[D.Angle]) -> Tensor[D.Dimless]:
    return tensor._base_tensor_class[D.Dimless]( # type: ignore
        np.cos(tensor._values)
    )

@overload
def sin(tensor: Scalar[D.Angle]) -> Scalar[D.Dimless]: ... # type: ignore

@overload
def sin(tensor: ScalarArray[D.Angle]) -> ScalarArray[D.Dimless]: ... # type: ignore

def sin(tensor: Tensor[D.Angle]) -> Tensor[D.Dimless]:
    return tensor._base_tensor_class[D.Dimless]( # type: ignore
        np.sin(tensor._values)
    )

@overload
def tan(tensor: Scalar[D.Angle]) -> Scalar[D.Dimless]: ... # type: ignore

@overload
def tan(tensor: ScalarArray[D.Angle]) -> ScalarArray[D.Dimless]: ... # type: ignore

def tan(tensor: Tensor[D.Angle]) -> Tensor[D.Dimless]:
    return tensor._base_tensor_class[D.Dimless]( # type: ignore
        np.tan(tensor._values)
    )

@overload
def acos(tensor: Scalar[D.Dimless]) -> Scalar[D.Angle]: ... # type: ignore

@overload
def acos(tensor: ScalarArray[D.Dimless]) -> ScalarArray[D.Angle]: ... # type: ignore

def acos(tensor: Tensor[D.Dimless]) -> Tensor[D.Angle]:
    return tensor._base_tensor_class[D.Angle]( # type: ignore
        np.acos(tensor._values)
    )

@overload
def asin(tensor: Scalar[D.Dimless]) -> Scalar[D.Angle]: ... # type: ignore

@overload
def asin(tensor: ScalarArray[D.Dimless]) -> ScalarArray[D.Angle]: ... # type: ignore

def asin(tensor: Tensor[D.Dimless]) -> Tensor[D.Angle]:
    return tensor._base_tensor_class[D.Angle]( # type: ignore
        np.asin(tensor._values)
    )

@overload
def atan(tensor: Scalar[D.Dimless]) -> Scalar[D.Angle]: ... # type: ignore

@overload
def atan(tensor: ScalarArray[D.Dimless]) -> ScalarArray[D.Angle]: ... # type: ignore

def atan(tensor: Tensor[D.Dimless]) -> Tensor[D.Angle]:
    return tensor._base_tensor_class[D.Angle]( # type: ignore
        np.atan(tensor._values)
    )

@overload
def atan2(y: ScalarArray[SomeDim], x: ScalarArray[SomeDim]) -> ScalarArray[D.Angle]: ...

@overload
def atan2(y: Scalar[SomeDim], x: Scalar[SomeDim]) -> Scalar[D.Angle]: ...

def atan2(y: Scalar | ScalarArray, x: Scalar | ScalarArray) -> Scalar | ScalarArray:
    return np.atan2(y.base_unit_value, x.base_unit_value) * Quantity.rad

def normalize_angle(angle: Scalar[D.Angle]) -> Scalar[D.Angle]:
    """Returns the given angle in its [0, 2π] range."""
    return angle % FULL_REV

def normalize_angle_symmetric(angle: Scalar[D.Angle]) -> Scalar[D.Angle]:
    """Returns the given angle in its [-π, π] range."""
    normalized = normalize_angle(angle)
    return normalized if normalized <= HALF_REV else normalized - FULL_REV

def angle2dms(angle: Scalar[D.Angle]) -> str:
    """Representation of the angle in DSM notation (degrees, minutes, seconds)

    Example: `39° 17′ N, 76° 36′ O`"""
    
    angle2convert = abs(angle.magnitude("deg"))
    deg, deg_dec = divmod(angle2convert, 1)
    min, min_dec = divmod(deg_dec * 60, 1)
    sec, _ = divmod(min_dec * 60, 1)
    
    return f"{deg}° {min}′ {sec}″"