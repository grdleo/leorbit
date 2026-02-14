
from abc import ABC, abstractmethod
from typing import Any, Generic, Self, TypeVar, cast

import numpy as np

from leorbit2.mathematics.matrix33 import Matrix33
from leorbit2.mathematics.dimensions import D, Dim, D, Number
from leorbit2.mathematics.scalar import Scalar
from leorbit2.mathematics.vector3 import Vector3


T1 = TypeVar("T1")
T2 = TypeVar("T2")
SomeDim = TypeVar("SomeDim", bound=Dim)

class Transform(Generic[T1, T2], ABC):
    """Abstract reversible transform between two tensor-like types."""

    @abstractmethod
    def do(self, tensor: T1) -> T2: ...

    @abstractmethod
    def undo(self, tensor: T2) -> T1: ...

    @abstractmethod
    def copy(self) -> Self: ...

    def reverse(self) -> Transform[T2, T1]:
        """Return a transform with ``do`` and ``undo`` swapped."""
        t = self.copy()

        do = t.do
        undo = t.undo

        tt = cast(Transform[T2, T1], t)

        tt.do = undo  # type: ignore
        tt.undo = do  # type: ignore

        return tt
    
class TransformIdentify(Generic[T1], Transform[T1, T1]):
    """Identity transform that leaves values unchanged."""

    def do(self, tensor: T1) -> T1:
        """Return the input value unchanged."""
        return tensor
    
    def undo(self, tensor: T1) -> T1:
        """Return the input value unchanged."""
        return tensor
    
    def copy(self) -> Self:
        """Return a new identity transform instance."""
        t = TransformIdentify[T1]()
        return cast(Self, t)
    
class TransformVector3Linear(Generic[SomeDim], Transform[Vector3[SomeDim], Vector3[SomeDim]]):
    """Linear transform for vectors using a dimensionless 3×3 matrix."""

    def __init__(self, matrix: Matrix33[D.Dimless]):
        """Initialize with the transformation matrix."""
        self.matrix = matrix
    
    def do(self, tensor: Vector3[SomeDim]) -> Vector3[SomeDim]:
        """Apply the linear transformation."""
        return self.matrix @ tensor
    
    def undo(self, tensor: Vector3[SomeDim]) -> Vector3[SomeDim]:
        """Apply the inverse linear transformation."""
        return self.matrix.inverse() @ tensor
    
    def copy(self) -> Self:
        """Return a deep copy of this transform."""
        t = TransformVector3Linear[SomeDim](
            matrix=self.matrix.copy()
        )
        return cast(Self, t)

class TransformVector3Affine(Generic[SomeDim], Transform[Vector3[SomeDim], Vector3[SomeDim]]):
    """Affine transform combining linear map and translation."""

    def __init__(self, matrix: Matrix33[D.Dimless], translation: Vector3[SomeDim]):
        """Initialize with matrix and translation components."""
        self.matrix = matrix
        self.translation = translation
    
    def do(self, tensor: Vector3[SomeDim]) -> Vector3[SomeDim]:
        """Apply affine transform ``M @ v + t``."""
        return self.matrix @ tensor + self.translation
    
    def undo(self, tensor: Vector3[SomeDim]) -> Vector3[SomeDim]:
        """Apply inverse affine transform ``M⁻¹ @ (v - t)``."""
        return self.matrix.inverse() @ (tensor - self.translation)
    
    def copy(self) -> Self:
        """Return a deep copy of this transform."""
        t = TransformVector3Affine[SomeDim](
            matrix=self.matrix.copy(), 
            translation=self.translation.copy()
        )
        return cast(Self, t)
    
class TransformVector3RotationZ(Generic[SomeDim], TransformVector3Linear[SomeDim]):
    """Specialized linear transform: rotation around the Z axis."""

    def __init__(self, angle_rad: Number | Scalar[D.Angle]):
        """Build the Z-rotation matrix from an angle in radians or angle scalar."""
        if isinstance(angle_rad, Scalar):
            angle_rad = angle_rad.magnitude("rad")
        
        c, s = np.cos(angle_rad), np.sin(angle_rad)
        rot_mat = Matrix33[D.Dimless].new(
            c, -s, 0,
            s, c, 0,
            0, 0, 1
        )

        super().__init__(rot_mat)
    
class TransformChain(Generic[T1, T2], Transform[T1, T2]):
    def __init__(self, *transforms: Transform[Any, Any]):
        """Create a composite transform executed in the given order."""
        self.transforms = list(transforms)

    def do(self, tensor: T1) -> T2:
        """Apply all transforms in forward order."""
        result: Any = tensor
        for t in self.transforms:
            result = t.do(result)
        return cast(T2, result)

    def undo(self, tensor: T2) -> T1:
        """Apply all inverse transforms in reverse order."""
        result: Any = tensor
        for t in reversed(self.transforms):
            result = t.undo(result)
        return cast(T1, result)
    
    def copy(self) -> Self:
        """Return a deep copy of the transform chain."""
        t = TransformChain[T1, T2](
            *(t.copy() for t in self.transforms)
        )

        return cast(Self, t)