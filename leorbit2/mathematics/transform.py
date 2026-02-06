
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
    @abstractmethod
    def do(self, tensor: T1) -> T2: ...

    @abstractmethod
    def undo(self, tensor: T2) -> T1: ...

    @abstractmethod
    def copy(self) -> Self: ...

    def reverse(self) -> Transform[T2, T1]:
        t = self.copy()

        do = t.do
        undo = t.undo

        tt = cast(Transform[T2, T1], t)

        tt.do = undo  # type: ignore
        tt.undo = do  # type: ignore

        return tt
    
class TransformIdentify(Generic[T1], Transform[T1, T1]):
    def do(self, tensor: T1) -> T1:
        return tensor
    
    def undo(self, tensor: T1) -> T1:
        return tensor
    
    def copy(self) -> Self:
        t = TransformIdentify[T1]()
        return cast(Self, t)
    
class TransformVector3Linear(Generic[SomeDim], Transform[Vector3[SomeDim], Vector3[SomeDim]]):
    def __init__(self, matrix: Matrix33[D.Dimless]):
        self.matrix = matrix
    
    def do(self, tensor: Vector3[SomeDim]) -> Vector3[SomeDim]:
        return self.matrix @ tensor
    
    def undo(self, tensor: Vector3[SomeDim]) -> Vector3[SomeDim]:
        return self.matrix.inverse() @ tensor
    
    def copy(self) -> Self:
        t = TransformVector3Linear[SomeDim](
            matrix=self.matrix.copy()
        )
        return cast(Self, t)

class TransformVector3Affine(Generic[SomeDim], Transform[Vector3[SomeDim], Vector3[SomeDim]]):
    def __init__(self, matrix: Matrix33[D.Dimless], translation: Vector3[SomeDim]):
        self.matrix = matrix
        self.translation = translation
    
    def do(self, tensor: Vector3[SomeDim]) -> Vector3[SomeDim]:
        return self.matrix @ tensor + self.translation
    
    def undo(self, tensor: Vector3[SomeDim]) -> Vector3[SomeDim]:
        return self.matrix.inverse() @ (tensor - self.translation)
    
    def copy(self) -> Self:
        t = TransformVector3Affine[SomeDim](
            matrix=self.matrix.copy(), 
            translation=self.translation.copy()
        )
        return cast(Self, t)
    
class TransformVector3RotationZ(Generic[SomeDim], TransformVector3Linear[SomeDim]):
    def __init__(self, angle_rad: Number | Scalar[D.Angle]):
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
        """executed in the given order"""
        self.transforms = list(transforms)

    def do(self, tensor: T1) -> T2:
        result: Any = tensor
        for t in self.transforms:
            result = t.do(result)
        return cast(T2, result)

    def undo(self, tensor: T2) -> T1:
        result: Any = tensor
        for t in reversed(self.transforms):
            result = t.undo(result)
        return cast(T1, result)
    
    def copy(self) -> Self:
        t = TransformChain[T1, T2](
            *(t.copy() for t in self.transforms)
        )

        return cast(Self, t)